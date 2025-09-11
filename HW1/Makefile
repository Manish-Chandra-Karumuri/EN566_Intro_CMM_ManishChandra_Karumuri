FXN=cos
WRITE_FILE=output.txt
READ_FILE=input.txt
FMT=png,jpeg,pdf,eps

.PHONY: plot print_plot write read print_read clean

plot:
	python trigonometry.py --function=$(FXN)

print_plot:
	@for ext in $$(echo $(FMT) | tr ',' ' '); do \
		python trigonometry.py --function=$(FXN) --print=$$ext; \
	done
	@echo "Plots for Functions saved!"

write:
	python trigonometry.py --function=$(FXN) --write=$(WRITE_FILE)

read:
	python trigonometry.py --read_from_file=$(READ_FILE)
	@echo "File read successful. Plot generated"

print_read:
	@for ext in $$(echo $(FMT) | tr ',' ' '); do \
		python trigonometry.py --read_from_file=$(READ_FILE) --print_txt_plt=$$ext; \
	done
	@echo "Plots from read file is saved!"

clean:
	@for ext in $$(echo $(FMT) | tr ',' ' '); do \
		rm -f Plots_from_Fxns.$$ext Data_Generated_Plot_From_Txt_file.$$ext; \
	done
	@rm -f $(WRITE_FILE) $(READ_FILE)
	@echo "Removed generated plots and text files."