docker run \
--rm \
-e JUPYTER_ENABLE_LAB=yes \
-v $PWD:/home/jovyan \
-p 8888:8888 \
geertvangeest/ngs-longreads-jupyter:latest \
start-notebook.sh

docker run \
--rm \
-v $PWD:/home/jovyan \
-p 8787:8787 \
-e PASSWORD=test \
geertvangeest/ngs-longreads-rstudio:latest