# acceleration_phewas
Scripts for UK Biobank accelerometer PheWAS analysis

# citation
See Khurshid et al. Accelerometer-Derived and Genetically Inferred Physical Activity and Human Disease. medRxiv 2021. https://doi.org/10.1101/2021.08.05.21261586

# overview
- scripts arranged into groups based on functionality
  - phecodes: scripts to generate phecode/outcome tables for association testing
  - seed_files: scripts to generate exposure tables for association testing
  - analysis/analysis: scripts to perform PheWAS analysis (relies on pre-built UK Biobank disease tables, see 'phecodes'), designed to run on outputs of 'phecodes' and 'seed_files'
  - analysis/plotting: scripts to process and visualize association tests
- 'covar' denotes fully adjusted models
