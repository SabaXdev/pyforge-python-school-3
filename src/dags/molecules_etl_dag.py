import boto3
from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine


# Default arguments for the DAG
default_args = {
    'owner': 'Saba',
    'depends_on_past': False,
    'email': ['sabakatamadze@gmail.com'],
    'email_on_failure': False,
    'email_on_retry': False,
    'start_date': days_ago(1),
    'retries': 1,
}

# Database connection setup
SQLALCHEMY_DATABASE_URL = "sqlite:///./molecules.db"

engine = create_engine(
    SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False}
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False,
                            bind=engine)


# 1. Extract Task
def extract_data(**kwargs):
    query = """
    SELECT name, description 
    FROM molecules 
    WHERE date(created_at) = date('now')
    """
    df = pd.read_sql(query, engine)
    return df.to_dict()


# 2. Transform Task
def transform_data(**kwargs):
    ti = kwargs['ti']
    data = ti.xcom_pull(task_ids='extract_data')
    df = pd.DataFrame(data)

    # RDKit descriptors
    def calculate_properties(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            h_donors = Descriptors.NumHDonors(mol)
            h_acceptors = Descriptors.NumHAcceptors(mol)
            lipinski_pass = (mw < 500 and logp < 5 and h_donors <= 5
                             and h_acceptors <= 10)
            return mw, logp, tpsa, h_donors, h_acceptors, lipinski_pass
        return None, None, None, None, None, False

    df[['Molecular Weight', 'LogP', 'TPSA', 'H Donors',
        'H Acceptors', 'Lipinski Pass']] = df['smiles'].apply(
        lambda s: pd.Series(calculate_properties(s))
    )
    return df.to_dict()


# 3. Load Task
def load_data_to_s3(**kwargs):
    """
    Save the transformed data to an Excel file and upload it to S3.
    """
    # First, pull the transformed data from XCom
    ti = kwargs['ti']
    transformed_data = ti.xcom_pull(task_ids='transform_data')

    # Convert the data into a DataFrame
    df = pd.DataFrame(transformed_data)

    # Save DataFrame as .xlsx
    file_path = '/tmp/molecule_data.xlsx'
    df.to_excel(file_path, index=False)

    # Upload the file to AWS S3 using Boto3
    s3_client = boto3.client('s3')
    bucket_name = 'hw-bucket-8844'  # Replace with your actual S3 bucket name
    object_name = 'molecule_data.xlsx'  # Name of the file in S3

    # Upload the file to S3 bucket
    s3_client.upload_file(file_path, bucket_name, object_name)


# Define the DAG
with DAG(
        'molecule_data_pipeline',
        default_args=default_args,
        description='Extract, Transform, Load molecules data to S3',
        schedule_interval='@daily',
        catchup=False,
) as dag:

    # Task Definitions
    extract_task = PythonOperator(
        task_id='extract_data',
        python_callable=extract_data
    )

    transform_task = PythonOperator(
        task_id='transform_data',
        python_callable=transform_data
    )

    load_task = PythonOperator(
        task_id='load_data',
        python_callable=load_data_to_s3
    )

    # Set Task Dependencies
    extract_task >> transform_task >> load_task
