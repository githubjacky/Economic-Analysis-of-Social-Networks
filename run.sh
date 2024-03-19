#!/bin/sh

docker run --rm \
	--name jupyter \
	--network host \
	--cap-add CAP_SYS_PTRACE \
	--shm-size 20gb \
	-e DOCKER_USER=$DOCKER_USER \
	-e UCX_TLS=tcp,cuda_copy \
	-e GRAPHISTRY_API_KEY=$GRAPHISTRY_API_KEY \
	-v $PWD/course:/home/$DOCKER_USER/economic-analysis-of-social-networks/course \
	-v $PWD/hw:/home/$DOCKER_USER/economic-analysis-of-social-networks/hw \
	--gpus '"device=3,capabilities=gpu"' \
	-p 8888:8888 \
	0jacky/economic-analysis-of-social-networks:latest \
	jupyter lab --ip 0.0.0.0 --port 8891 --no-browser --allow-root
