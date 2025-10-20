PYTHON = python
SCRIPTS = oscillator.py poisson.py rnumbers.py rwalk.py diffusion.py gases.py
PART ?= 1,2,3,4,5

all: oscillator poisson rnumbers rwalk diffusion gases

oscillator: oscillator.py
	@echo "--- Running Oscillator Simulation (Part(s): $(PART)) ---"
	$(PYTHON) oscillator.py --part=$(PART)

poisson: poisson.py
	@echo "--- Running Poisson Equation Solver (Part(s): $(PART)) ---"
	$(PYTHON) poisson.py --part=$(PART)

rnumbers: rnumbers.py
	@echo "--- Running Random Number Analysis (Part(s): $(PART)) ---"
	$(PYTHON) rnumbers.py --part=$(PART)

rwalk: rwalk.py
	@echo "--- Running Random Walk Simulation (Part(s): $(PART)) ---"
	$(PYTHON) rwalk.py --part=$(PART)

diffusion: diffusion.py
	@echo "--- Running Diffusion Equation Solver ---"
	$(PYTHON) diffusion.py --part=2

gases: gases.py
	@echo "--- Running Gas Mixing Simulation (Part(s): $(PART)) ---"
	$(PYTHON) gases.py --part=$(PART)

clean:
	@echo "Cleaning up generated plot and animation files..."
	rm -f *.png *.gif
	@echo "Cleanup complete."

.PHONY: all clean oscillator poisson rnumbers rwalk diffusion gases

