This folder contains all data processing and analysis scripts 

## Main Files

### Data Cleaning
- **clean.do**: Primary data cleaning script that processes raw whaling voyage data, creates the master dataset, and merges captain information with voyage records. This must be run first.

### Initial Exploration Files
- **IE[1-9]*.do/r**: Sequential analysis files exploring different aspects of the data. These should be run in numerical order after clean.do

## Subfolders

### ML_estimation/
Contains Maximum Likelihood estimation code for structural production models with captain and vessel fixed effects. **This folder has its own README with detailed documentation of the estimation procedures.**

## Output

All code files generate outputs following these conventions:
- Figures are saved to `../Figures/` with prefix matching the source file
- Tables are saved to `../Writeup/Tables/` in LaTeX format
- Intermediate data files (temp*) are saved to `../Data/`