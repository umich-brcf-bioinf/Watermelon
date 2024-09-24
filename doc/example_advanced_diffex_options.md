# Example - Advanced Diffex Options

In the README, a straightforward differential expression setup was covered. This document outlines some additional possibilities of configuration that are helpful for more advanced differential expression analysis setups.

> Note: For maximum flexibility, the majority of these advanced options have only basic error checking built around them. It's best to double-check all of your configuration details against the documentation and/or code, and ensure that your combination of configuration options makes sense.

## Using `subset`:

It is possible to subset the samples based on a samplesheet characteristic while performing the differential expression. The samples will be included in plots upstream of differential expression analysis e.g. expression QC plots, but will be excluded before creating the DESeqDataset. This is useful when you desire this combination. If you want to exclude samples from the expression QC plots *and* the differential expression analysis, then removing them from the samplesheet (excluding them entirely) is recommended.

An example specifying `subset` in the config, so that only samples with `high_fat` in the `diet` column of the samplesheet are included during differential expression analysis, looks like this:

```
diffex:
    count_min_cutoff: 1 #Filter lowly expressed genes using this cutoff
    model_genotype_subset_high_fat:
        subset: diet::high_fat
        linear_fold_change: 1.5
        adjustedPValue: 0.05
        contrasts:
        -   ko_v_wt
        DESeq2:
            design: ~ genotype
            factor_name: genotype
            DESeq:
                betaPrior: true
            results:
                contrast: list(paste0(factor_name, test_name), paste0(factor_name, reference_name))
```

## Using `set_ref_levels`:

The typical way of defining contrasts as outlined in the README walkthrough is conveyed explicitly in the config by the contrast string, e.g. `test_v_reference`. However, sometimes we want to perform analyses which specify things in a different way - the LRT example below is one instance. For a given attribute, DESeq will use the first item, sorted alphabetically, by default as the reference level when pulling results. In many cases this is not what we want, so in the config it is possible to explicitly specify a reference level for each column in the samplesheet. Provide a list of `factor::level` pairs and each in the list will set the reference level of `factor` as `level`, e.g.:

```
set_ref_levels:
  - treatment::control
  - diet::normal
```

## LRT example:

One common type of analysis specifies the `LRT` test for `DESeq2::DESeq()`. Here is an example config section for that, with some bullet points below highlighting some useful patterns that may not be obvious at first glance:

```
diffex:
    count_min_cutoff: 1
    model_treatment:
        adjustedPValue: 0.05
        contrasts:  # Empty contrasts triggers different behavior
        set_ref_levels:
          - treatment::control
        DESeq2:
            design: ~ batch + treatment
            DESeq:
                test: LRT
                reduced: ~ treatment
            results: # Leave blank to use default args for results
```

- Notice `linear_fold_change` is not specified.
  - By excluding it or leaving it blank, fold-change will not be considered when calling differentially expressed genes, only P-val
- Notice `contrasts` is left blank
  - By excluding it or leaving it blank, a different behavior is triggered vs when it is populated
  - Runs through DESeq() and results() just once, creating output files with default names like `results.annot.txt`
  - Does not produce volcano plots
  - Arguments to `results` section not dependent on `contrasts` or `factor_name`
- Notice `factor_name` is not specified. It is not needed, since we are not using contrast strings or `factor_name` to specify arguments to `results`.
- Notice `set_ref_levels` is used. This functionality is explained in more detail in its own section of this document
- Notice `results` is left blank. Take note that in this case, an empty `results` section has different behavior than excluding `results` entirely
  - When the `results` section is left blank, it will use default arguments for the call to `DESeq2::results`
  - Of course, you always have the option of specifying arguments for `results` as needed
  - This is just a simple shortcut for the common use case of using default values here
  - Excluding the `results()` section entirely (incl. the `results` key itself) should only be used in conjunction with `lfcShrink()` as described below

## Using `lfcShrink()` instead of `results()`:

It is possible to specify parameters for `DESeq2::lfcShrink()` and call that function instead of `DESeq2::results()` to pull out results. This section and the `results` section are mutually exclusive. An example config is here:

```
diffex:
    count_min_cutoff: 1 #Filter lowly expressed genes using this cutoff
    model_treatment_inWT:
        linear_fold_change: 1.5
        adjustedPValue: 0.05
        contrasts:
        -   RM1_v_noRM1
        DESeq2:
            design: ~ Mouse + Treatment
            factor_name: Treatment
            DESeq:
                betaPrior: false
            #The below parameters are also parsed by R and used to call DESeq2::lfcShrink()
            lfcShrink:
                 coef: paste(factor_name, test_name, 'vs', reference_name, sep='_')
                 type: apeglm
```
