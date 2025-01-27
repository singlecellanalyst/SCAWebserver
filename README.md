Dockers for the SingleCellAnalyst webserver: www.singlecellanalyst.org

Instructions for Using the Docker Version of the SingleCellAnalyst:

1. Prerequisites:
- Ensure that Docker is installed on your system. If not, download and install Docker from https://www.docker.com/.
- Verify Docker installation by running the command on the terminal:
```sh
docker --version
```
2. Pull the Docker Image:
- Download the Docker image of your selected omics analysis pipeline from the repository:
```sh
https://github.com/singlecellanalyst?tab=packages&repo_name=SCAWebserver
```

3. Run the Docker Container:
- Launch the platform by executing:
```sh
docker run -p 8080:8080 -d <the docker name of the omics selected>
```

4. Access the Platform:
- Open a web browser and go to:
```sh
http://localhost:8080
```
The SingleCellAnalyst interface of the omics type you have downloaded should load, allowing you to upload datasets and start your analyses.

5. Stop the Docker Container (Optional):
- To stop the running Docker container, execute:
```sh
docker stop <container id>
```
- You can find the <container_id> by running:
```sh
docker ps
```
