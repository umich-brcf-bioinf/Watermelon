## Git Hooks for Automated Tests

For this repository, we use git hooks to run the appropriate automated tests locally. This provides a suitable solution for testing bioinformatics pipelines, which need large files and relatively long processing times in order to adequately test. This solution allows us to perform complete test-runs of the pipeline on the example dataset, using the cluster to perform jobs in an identical way to what we use on a day-to-day basis.

Since git hooks are not preserved in git repositories, they must be manually added to each cloned repository. There are two files added to our `.git/hooks/` directory to automatically invoke this testing functionality. They have identical contents, set up to handle direct commits or merges to production in the same fashion.

Contents of `pre-commit` & `pre-merge-commit` files are identical:

    #! /bin/sh

    BRANCH=$(git rev-parse --abbrev-ref HEAD)

    if [ $BRANCH == "master" ]; then
        if [[ $HOSTNAME =~ "gl-login" ]]; then
            pytest --basetemp /scratch/cgates_root/cgates1/trsaari/tmp/ tests/
        else
            echo "Must be on HPC to run functional tests" &&
            exit 1
        fi
    else
        pytest --ignore=tests/snakemake_functional_run_test.py tests/
    fi


This means that the lightweight unit tests will be performed whenever a commit is made e.g. on a development branch, so regular development can happen anywhere. However, whenever something is commited or merged into production, we verify that we're on one of the HPC login nodes and then run the full suite of tests.

**Note: You will have to set `--basetemp` to a directory that you have access to**

**Note: The `gl-login` string in the $HOSTNAME variable works for the Great Lakes HPC at UMich. This may need to be updated as well**
