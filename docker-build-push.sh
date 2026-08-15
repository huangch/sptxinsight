docker build -f ./Dockerfile -t sptxinsight:latest .
docker tag sptxinsight:latest huangchtw/sptxinsight:latest
docker push huangchtw/sptxinsight:latest


