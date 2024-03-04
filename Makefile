.PHONY: jupyter


jupyter:
	docker compose run --rm --service-ports jupyter-lab
