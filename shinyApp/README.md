# Shiny app

This is the code for a [Shiny](https://shiny.posit.co/) app that allows you to
explore the data in this study interactively.

## Docker image

A docker image of this app you can run locally can be found [here](https://hub.docker.com/repository/docker/multitude5286/ad-bbb-shinyapp/general), steps to use:

1. Install Docker, documentation [here](https://docs.docker.com/desktop/)
2. Pull the container
3. Run the container, and view in your web browser at `http://localhost:3838`

```bash
# pull the container
docker pull multitude5286/ad-bbb-shinyapp:latest
# run the container
docker run --rm -p 3838:3838 multitude5286/ad-bbb-shinyapp
```
