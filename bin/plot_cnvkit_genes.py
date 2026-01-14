#!/usr/bin/env python3

import pandas as pd
import gzip
import matplotlib.pyplot as plt
import argparse


def main(args):
    gene = args.gene

    gtf_file = args.gtf
    cns_file = args.cns
    cnr_file = args.cnr

    out_file = f"{args.gene}.png"
    if args.output:
        out_file = args.output

    plot_gene_cnv(
        gene,
        gtf_file,
        cnr_file,
        cns_file,
        flank=50_000,
        ymin=-5,
        ymax=5,
        out_png=out_file,
        vcf_file=args.scored_vcf,
    )


def plot_gene_cnv(
    gene,
    gtf_file,
    cnr_file,
    cns_file,
    flank=50_000,
    ymin=-5,
    ymax=5,
    out_png=None,
    vcf_file=None,
):
    # Load GTF
    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "chrom",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attr",
        ],
        low_memory=False,
    )

    # Only get data for supplied gene
    gtf["gene_name"] = gtf["attr"].str.extract(r'gene_name "([^"]+)"')

    gene_gtf = gtf[(gtf["gene_name"] == gene)]

    if gene_gtf.empty:
        raise ValueError(f"Gene '{gene}' not found in GTF")

    chrom = gene_gtf.iloc[0]["chrom"]

    # make sure is string, for normalisation of chr
    chrom = str(chrom)

    # GTF is 1-based cnvkit bedlike is not, convert
    gene_start = gene_gtf["start"].min() - 1
    gene_end = gene_gtf["end"].max()

    region_start = max(0, gene_start - flank)
    region_end = gene_end + flank

    exons = gene_gtf[gene_gtf["feature"] == "exon"]

    # Load CNVkit data
    cnr = pd.read_csv(cnr_file, sep="\t", low_memory=False)
    cns = pd.read_csv(cns_file, sep="\t", low_memory=False)

    # make sure it is string
    cnr["chromosome"] = cnr["chromosome"].astype(str)
    cns["chromosome"] = cns["chromosome"].astype(str)

    # Normalize chr prefix based on CNVkit files
    cnv_has_chr = cnr["chromosome"].str.startswith("chr").any()

    if cnv_has_chr and not chrom.startswith("chr"):
        chrom = "chr" + chrom
    elif not cnv_has_chr and chrom.startswith("chr"):
        chrom = chrom.replace("chr", "")

    cnr_r = cnr[
        (cnr["chromosome"] == chrom)
        & (cnr["start"] < region_end)
        & (cnr["end"] > region_start)
    ]

    cns_r = cns[
        (cns["chromosome"] == chrom)
        & (cns["start"] < region_end)
        & (cns["end"] > region_start)
    ]

    # Plot segments, probes and exons
    fig, ax = plt.subplots(figsize=(10, 4))

    # CNV probes/bins
    ax.scatter(
        (cnr_r["start"] + cnr_r["end"]) / 2,
        cnr_r["log2"],
        s=10,
        alpha=0.5,
        color="grey",
    )

    # called segments
    for _, seg in cns_r.iterrows():
        ax.hlines(seg["log2"], seg["start"], seg["end"], linewidth=2, color="firebrick")

    # exons
    exon_y = ymin + 0.2
    exon_h = 0.5

    for _, exon in exons.iterrows():
        ax.add_patch(
            plt.Rectangle(
                (exon["start"] - 1, exon_y),
                exon["end"] - exon["start"],
                exon_h,
                color="steelblue",
            )
        )

    # introns in between
    exons_sorted = exons.sort_values("start")

    for i in range(len(exons_sorted) - 1):
        x1 = exons_sorted.iloc[i]["end"]
        x2 = exons_sorted.iloc[i + 1]["start"] - 1

        ax.hlines(
            y=exon_y + exon_h / 2, xmin=x1, xmax=x2, color="steelblue", linewidth=1
        )

    # highlight gene
    ax.axvspan(gene_start, gene_end, color="steelblue", alpha=0.1)

    sv_ymin = ymin  # bottom of CNV plot
    sv_ymax = exon_y + exon_h + 0.2  # just below exons

    if vcf_file:
        sv_df = parse_vcf(vcf_file, chrom, region_start, region_end)
        
        if not sv_df.empty:
            sv_df = sv_df.sort_values("start")

            tracks = []
            text_base_y = sv_ymax + 0.2
            text_step = 0.35 

            for _, sv in sv_df.iterrows():
                ax.axvspan(
                    sv["start"],
                    sv["end"],
                    ymin=(sv_ymin - ymin) / (ymax - ymin),
                    ymax=(sv_ymax - ymin) / (ymax - ymin),
                    color="orange",
                    alpha=0.3,
                )

                # try to keep track of variants overlapping, avoid text jumble
                for track_idx, last_end in enumerate(tracks):
                    if sv["start"] > last_end:
                        tracks[track_idx] = sv["end"]
                        break
                else:
                    track_idx = len(tracks)
                    tracks.append(sv["end"])

                ax.text(
                    (sv["start"] + sv["end"]) / 2,
                    text_base_y + track_idx * text_step,
                    f'{sv["SVTYPE"]}:{sv["RankScore"]}',
                    fontsize=7,
                    ha="center",
                    va="top",
                )

    ax.set_xlim(region_start, region_end)
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel("log2 copy ratio")
    ax.set_xlabel(f"{chrom}:{region_start:,}-{region_end:,}")
    ax.set_title(gene)

    plt.tight_layout()

    if out_png:
        plt.savefig(out_png, dpi=200)
        plt.close(fig)
    else:
        plt.show()


def parse_vcf(vcf_file, chrom, region_start, region_end):
    """
    Returns a dataframe of SV calls overlapping the given region.
    Expects INFO field with SVTYPE and RankScore.
    """
    svs = []
    open_func = gzip.open if vcf_file.endswith(".gz") else open

    with open_func(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            v_chrom, pos, _, _, _, _, _, info = fields[:8]
            if str(v_chrom) != str(chrom):
                continue
            pos = int(pos)

            # parse INFO field
            info_dict = dict(item.split("=") for item in info.split(";") if "=" in item)
            svtype = info_dict.get("SVTYPE", "NA")
            rankscore = info_dict.get("RankScore", "NA")

            # end is not VCF standard but we use it, fallback
            end = int(info_dict.get("END", pos))
            
            if not (pos <= region_end and end >= region_start):
                continue

            svs.append(
                {"start": pos, "end": end, "SVTYPE": svtype, "RankScore": rankscore}
            )

    return pd.DataFrame(svs)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plots gene coverage for provided gene name for a given CNVkit called sample"
    )

    parser.add_argument(
        "--cnr",
        "-r",
        type=str,
        required=True,
        help="Path to CNVkit regions file (.cnr)",
    )
    parser.add_argument(
        "--cns",
        "-s",
        type=str,
        required=True,
        help="Path to CNVkit segments file (.cns)",
    )
    parser.add_argument(
        "--gtf",
        "-d",
        type=str,
        required=True,
        help="path to gtf file, preferably mane transcripts only",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="name of output png, if not provided output will be gene_name.png",
    )
    parser.add_argument(
        "--gene", "-g", type=str, required=True, help="HGNC symbol name"
    )
    parser.add_argument(
        "--scored_vcf",
        "-v",
        type=str,
        help="Optionally add scored vcf to highlight calls",
    )
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
