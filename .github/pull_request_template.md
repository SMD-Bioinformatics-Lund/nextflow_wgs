<!--
Thanks for contributing to the CMD nextflow_wgs pipeline!

Please use the checklists below to document changes and performed tests.

Remember that this template doubles as our test/review documentation. 
-->
# Description and reviewer info

What is changed? How does this update improve the pipeline? (For reviewers) How to test it?


## Type of change
<!--
    Major change counts as a change that breaks backward compatibility
    Minor change is a substantial change that requires testing before deployment
    Patch is a minor change like a bug fix, code comment/style fix, etc.
-->
- [ ] Documentation
- [ ] Patch
- [ ] Minor change
- [ ] Major change 

# Checklist

- [ ] Self-review of my code
- [ ] Update the CHANGELOG
- [ ] Tag the latest commit (vX.Y.Z format)

<!--
    Select a checklist below based on selection under # Type of change
    and delete the sections that do not apply to this PR:
-->

## Documentation
- [ ] At least one other person has reviewed my changes (not required for trivial changes)

## Patch
- [ ] Stub run completes without errors or new warnings
- [ ] At least one other person has reviewed and approved my code (not required for trivial changes)

## Major / Minor change
- [ ] Stub run completes without errors or new warnings
- [ ] GIAB single (wgs) finishes and differences to current master branch have been investigated (using [PipeEval](https://github.com/Clinical-Genomics-Lund/PipeEval))
- [ ] GIAB trio (wgs) finishes and differences to current master branch have been investigated
- [ ] Seracare (onco) sample finishes and differences to current master branch have been investigated
- [ ] Output from processes with altered code have been inspected (i.e. in the work folder: .command.sh, .command.log, .command.err and result files)
- [ ] All three samples can be loaded into Scout
- [ ] At least one other person has reviewed and approved my code
- [ ] I have made corresponding changes to the documentation (software versions, etc.)

# Test/review documentation

## Review performed by

- [ ] Alexander
- [ ] Jakob
- [ ] Paul
- [ ] Ryan
- [ ] Viktor

(Add if missing)

## Testing performed by

- [ ] Alexander
- [ ] Jakob
- [ ] Paul
- [ ] Ryan
- [ ] Viktor
