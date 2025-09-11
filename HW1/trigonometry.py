import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


parser = argparse.ArgumentParser(description='Takes input functions and operations and execute them for the user')
parser.add_argument('--function', help = 'Functions to evaluate', required = False)
parser.add_argument('--read_from_file', help = 'Read from the file for input and calculate the output', required=False)
parser.add_argument('--write', help = 'Writes the input and output to a txt file', required = False)
parser.add_argument('--print', help = 'Saves the plots generated from user input', required = False)
parser.add_argument('--print_txt_plt', help = 'Saves the plot generated from the .txt file', required = False)

args = parser.parse_args()


x = np.arange(-10, 10 + 0.05, 0.05)


def func(x):
    if args.function:
        if "," in args.function:
            new_args = args.function.split(',')
            for i in range(len(new_args)):
                fnc = f"np.{new_args[i]}(x)"  
                y = eval(fnc)
                plt.plot(x, y, label = f'np.{new_args[i]}(x)')
            
            plt.xlabel('x values')
            plt.title(f'{new_args}(x) vs x values')
            plt.grid(True)
            plt.legend()
            plt.ylabel(f'{new_args}(x)')

        else:
            fn = f"np.{args.function}(x)"  
            y = eval(fn)
            plt.plot(x, y, label = f'{fn}')
            plt.legend()
            plt.xlabel('x values')
            plt.ylabel(f'{fn}')
            plt.title(f'{fn} vs x values')
            plt.grid(True)

        if args.print is not None:
            if ',' in args.print:
                formats = args.print.split(',')
                for fmt in formats:
                    plt.savefig(f'Plots_from_Fxns.{fmt}')
            else:
                plt.savefig(f'Plots_from_Fxns.{args.print}')
        
        plt.show()

def read_from(textfile):
    data = pd.read_csv(textfile, sep = '\t', header=None)
    
    
    if data.empty:
        print(f"Error: The file '{textfile}' is empty or could not be read properly.")
        return 

    x_values = data.columns[0]
    y_values = data.columns[1:]
    data[x_values] = pd.to_numeric(data[x_values], errors='coerce')
    x_val = data[0]
    x_first = x_val.iloc[1]
    x_last = x_val.iloc[-1]
    
    for y in y_values:
        data[y] = pd.to_numeric(data[y], errors='coerce')
    for y in y_values:
        plt.plot(data[x_values], data[y], label=f'Column data {y}')
    #print(x_first, x_last)
    plt.xlim(x_first, x_last)
    plt.grid(True)
    plt.title('Plots from the txt file input')
    plt.xlabel("x-values from file")
    plt.ylabel("y-values from file")
    plt.legend()

    if args.print_txt_plt is not None:
        if ',' in args.print_txt_plt:
            formats = args.print_txt_plt.split(',')
            for fmt in formats:
                plt.savefig(f'Data_Generated_Plot_From_Txt_file.{fmt}')
        else:
            plt.savefig(f'Data_Generated_Plot_From_Txt_file.{args.print_txt_plt}')

    plt.show()
    
def write_to(functions, textfile):
    data = pd.DataFrame({'x': x})
    funcs = [f.strip() for f in functions.split(',')]
    for func in funcs:    
        f = getattr(np, func)
        y = f(x)
        data[func] = y
    
    data.to_csv(textfile, sep = '\t', index=False)
    print(f"Data written to {textfile}")

if __name__== '__main__':
    if args.read_from_file:
        read_from(args.read_from_file)
    elif args.function:
        func(x)
        if args.write:
            write_to(args.function, args.write)
    else:
        print("No function or file specified. Plotting default function: 'cos'")
        args.function = 'cos'
        func(x)

