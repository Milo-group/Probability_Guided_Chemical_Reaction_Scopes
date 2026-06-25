This folder contains the following subfolders:  
  
**Data**  
All relevant ```*.csv``` files used for modeling, divided into three-class and binary folders.  
    
**Model reports**  
PDF files containing model-search reports for the 10% holdout classification (binary and 3-class ordinal):
- `Model.Search.Pyr.Phos.Binary.10pctHoldout.pdf`
- `Model.Search.Pyr.Phos.Ordinal.10pctHoldout.pdf`

Generate these quickly (without re-running `sub_model_log`) using `Generate_Model_Search_Reports.R`, which uses the pre-ranked model tables recorded in `McNally_10_percent_random.R` and the PDF helpers in `model_search_report_helpers.R`.

An ```*.R``` file is also included, allowing the recreation of the full analysis (including model search).
