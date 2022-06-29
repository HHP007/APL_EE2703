"""
NAME:HARIHARAN P
ROLL:EE20B042
EE2703_A1

"""
import sys # first import the required library sys
# Function to Extract the words/tokens from a Line
def linetowords(Line):
    total_Words = Line.split()

    #  if the file has only R, L, C, Independent Sources
    if(len(total_Words) == 4):
        elementName = total_Words[0]
        node1 = total_Words[1]
        node2 = total_Words[2]
        value = total_Words[3]
        return [elementName, node1, node2, value]

    # if the file has Current Controlled current/voltage Source
    elif(len(total_Words) == 5):
        elementName = total_Words[0]
        node1 = total_Words[1]
        node2 = total_Words[2]
        voltageSource = total_Words[3]
        value = total_Words[4]
        return [elementName, node1, node2, voltageSource, value]

    # if the file has Voltage Contolled current/voltage Source
    elif(len(total_Words) == 6):
        elementName = total_Words[0]
        node1 = total_Words[1]
        node2 = total_Words[2]
        voltageSourceNode1 = total_Words[3]
        voltageSourceNode2 = total_Words[4]
        value = total_Words[5]
        return [elementName, node1, node2, voltageSourceNode1, voltageSourceNode2, value]

    else:
        return [] #if no cases are satisfied , return an empty list

def printasreqorder(LinesTokens): #to print the tokens in reverse order
    for x in LinesTokens[::-1]:
        for y in x[::-1]:
            print(y, end=' ')
        print('')
    print('')
    return

from sys import argv, exit
# To improvise editablity
STARTING_CIR = ".circuit"
ENDING_CIR = ".end"
   
    # checking number of command line arguments
if len(sys.argv)!=2 :
        sys.exit("Invalid number of inputs! Pass the netlist file as the command line argument.")
else:
        try:
            actual_circuit_details = sys.argv[1]

            # checking if given netlist file is of correct type
            if (not actual_circuit_details.endswith(".netlist")):
                print("incorrect file type! please give .netlist file only")
            else:
                with open (actual_circuit_details, "r") as f:
                    SPICE_Lines = []
                    for line in f.readlines():
                        SPICE_Lines.append(line.split('#')[0].split('\n')[0])
                    try:
                        # finding the location of the identifiers
                        identifier1 = SPICE_Lines.index(STARTING_CIR)
                        identifier2 = SPICE_Lines.index(ENDING_CIR)

                        SPICELines_Actual = SPICE_Lines[identifier1+1:identifier2]
                        SPICELines_Tokens = [linetowords(line) for line in SPICELines_Actual]

                        # function to output Circuit details in Reverse Order
                        printasreqorder(SPICELines_Tokens)
                    except ValueError:
                        print("The Netlist given does not follow to the given format! Make sure to have .circuit and .end lines in the file.")
        except FileNotFoundError:
            print("Given file does not exist! Make sure you have entered the name of the netlist file correctly.")