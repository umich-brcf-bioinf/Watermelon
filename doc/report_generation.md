# Report Generation

Upon successful completion of the pipeline, two versions of a draft report will be generated in the analysis folder, as `report/report_draft.md` and `report/report_draft.html`. These are markdown and HTML files, respectively. The HTML version of the draft report can be opened and viewed in an internet browser, and the markdown version can be edited in a text editor.

The following steps describe the process of creating a final report:

1. Open the HTML report in an internet browser and inspect it for portions which may need editing. For example, background information may need to be manually entered or updated, certain figures or results may necessitate comments for context. Additionally, some figures or details may be inaccurate or redundant and thus should be removed.

2. With the above information in mind, edit the markdown version of the draft report, `report_draft.md`.

3. With the edited `report_draft.md` saved, generate a finalized report. To do this you'll call snakemake with the specified target `report_finalize`. E.g.:

    ```
    snakemake --configfile config_20190822.yaml --snakefile Watermelon/align_qc.smk \
        --profile Watermelon/config/profile-greatlakes report_finalize
    ```

    This will create the output `report/report_final.html`, and it will also copy this file and a corresponding multiQC report into `deliverables/report/`.

    Note: `multiqc_linked_report.html` should be kept co-located with the report, in order for the linked multiQC report to work.

4. Open `report_final.html` in an internet browser and inspect it. The edits that were just made to the markdown version should now be applied to the finalized HTML version. If anything did not turn out as planned, go back to `report_draft.md` and edit it as needed. Then re-run snakemake with the `report_finalize` target as in step 3.

    Note: If multiple rounds of editing are needed, it may be useful to store versions of `report_draft.md` with different filenames, so that the different versions can be compared or used together to create an ideal finalized report.

The final report, as noted, will reside in `deliverables/report/` as `report_final.html`, along with `multiqc_linked_report.html`
