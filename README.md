
# Molecules Project

## Overview

This project leverages GitHub Actions to automate the deployment of our application to an AWS EC2 instance whenever changes are pushed to the main or other branches. This ensures that our application is always up-to-date and running the latest code.

## GitHub Actions Workflow

The workflow is defined in a YAML file located in the `.github/workflows` directory. It contains steps for setting up the environment, installing dependencies, uploading code to the EC2 instance, and starting the application.

### 1. **Installing Dependencies**

When a push is made to the main branch, the workflow is triggered. The first step involves setting up the Python environment and installing necessary dependencies. This is done using the following steps:

```yaml
- name: Set up Python
  uses: actions/setup-python@v4
  with:
    python-version: '3.10'

- name: Install dependencies
  run: |
    python -m pip install --upgrade pip
    pip install -r src/requirements.txt
```

- **Set up Python:** Specifies the Python version required for the project.
- **Install dependencies:** Upgrades `pip` and installs all the required packages listed in the `requirements.txt` file.

### 2. **Uploading Code**

After the dependencies are installed, the workflow uploads the application code to the EC2 instance using SCP (Secure Copy Protocol). The relevant section of the workflow looks like this:

```yaml
- name: Upload code via SCP
  uses: appleboy/scp-action@master
  with:
    host: ${{ secrets.EC2_HOST }}
    username: ${{ secrets.EC2_USER }}
    key: ${{ secrets.EC2_SSH_KEY }}
    source: "./*"
    target: "/home/ubuntu/"
```

- **Host:** The public IP address or DNS of the EC2 instance, stored as a GitHub secret.
- **Username:** The username for the EC2 instance (e.g., `ubuntu`).
- **Key:** The SSH private key for accessing the EC2 instance, also stored as a GitHub secret.
- **Source:** The local files to upload, specified as `"./*"` to copy everything in the repository.
- **Target:** The directory on the EC2 instance where the code will be uploaded.

### 3. **Starting the Application**

Once the code is uploaded, the workflow SSHs into the EC2 instance to start the application. This is achieved through the following step:

```yaml
- name: SSH into EC2 instance and deploy
  uses: appleboy/ssh-action@master
  with:
    host: ${{ secrets.EC2_HOST }}
    username: ${{ secrets.EC2_USER }}
    key: ${{ secrets.EC2_SSH_KEY }}
    script: |
      sudo apt update && sudo apt install -y python3-pip
      pip3 install -r requirements.txt
      nohup python3 main.py &
```

- **Script Execution:** The script updates the package manager, installs `python3-pip`, installs the application requirements, and starts the application in the background using `nohup`.

### Conclusion

This GitHub Actions workflow streamlines the deployment process, ensuring that every push to the main branch results in an automatic deployment to our AWS EC2 instance. This not only saves time but also minimizes the chances of human error during the deployment process.
