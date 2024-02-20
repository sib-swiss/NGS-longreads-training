# got all data from PRJNA935757 delivered to SIB in the S3 bucket sra-data-delivery-sib

aws s3 sync --source-region us-east-1 s3://sra-data-delivery-sib .
