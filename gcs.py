# Import the required libraries, classes and modules
import pandas as pd
import numpy as np
from sdv.metadata import MultiTableMetadata
from sdv.single_table import GaussianCopulaSynthesizer

# Define the sample size
sample_size = 642

# Define the the number of sinthetic samples
n_samples = 10

# Read data from the CSV file into a DataFrame named 'db'
db = pd.read_csv("f.csv")

# Drop columns with too many levels from the DataFrame
db.drop(columns=['id', 'COUNTRY', 'SiteKey'], inplace=True)

# Assign a new 'id' column with values ranging from 0 to the length of the DataFrame
db = db.assign(id=range(len(db)))

# Create an instance of MultiTableMetadata
metadata = MultiTableMetadata()

# Detect metadata from the provided DataFrame 'db'
metadata.detect_from_dataframes(data={"db": db})

# Save the metadata to a JSON file named 'metadata.json'
metadata.save_to_json('metadata.json')

# Create an instance of GaussianCopulaSynthesizer
# The metadata.tables['db'] provides the metadata information for the 'db' table
# Verbose mode is enabled for additional output information
synthesizer = GaussianCopulaSynthesizer(metadata.tables['db'], verbose=True)

# Train the CTGAN synthesizer on the data in the 'db' DataFrame
synthesizer.fit(db)

# Import the evaluate_quality function from the sdv.evaluation.single_table module
from sdv.evaluation.single_table import evaluate_quality

# Evaluate the quality of the synthetic data
# 'db': Original data, synthesizer.sample(num_rows=sample_size): Synthetic data sample,
# metadata.tables['db']: Metadata information for the 'db' table
quality_report = evaluate_quality(
    db,
    synthesizer.sample(num_rows=sample_size),
    metadata.tables['db']
)

# Loop through n_samples iterations
for i in range(n_samples):
  
    # Generate synthetic data with sample_size rows
    synthetic_data_i = synthesizer.sample(num_rows=sample_size)
    
    # Define the file name with the iteration number
    file_name = "~/temp_{}.csv".format(i)
    
    # Save the synthetic data to a CSV file
    synthetic_data_i.to_csv(file_name, index=False)
