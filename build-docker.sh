docker build --build-arg UID=$(id -u) --build-arg GID=$(id -g) -f ./Dockerfile -t sptxinsight:latest .
docker tag sptxinsight:latest lj-docker-reg.pfizer.com/huangc78/sptxinsight
docker push lj-docker-reg.pfizer.com/huangc78/sptxinsight
