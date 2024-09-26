import boto3
import json

# Initialize my AWS credentials
session = boto3.Session(region_name='eu-north-1')
lambda_client = session.client('lambda')

# Define the event to send
event = {
    'names': ['Saba', 'Davit', 'Nikoloz', 'Gio']
}

# Invoke the Lambda function
response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event)
)

response_payload = json.loads(response['Payload'].read())
print(response_payload)
