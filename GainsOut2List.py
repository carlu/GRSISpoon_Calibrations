#!/usr/bin/python
import fileinput

for line in fileinput.input():
   Items = line.split()
   if (len(Items) == 11): 
      if(len(Items[0])==10):
         Name=Items[0]
      else:
         Name=Items[0]+"x"
      
      #print Items
      #print Name
      
      #string = Name + "   {\n"
      #string += "   NAME: " + Name + "\n"
      #string += "   ENG_COEFF: " + Items[2] + " " + Items[5] + " " + Items[8]
      #string += "\n}\n"
      
      string = Name + " " + Items[2] + " " + Items[5] + " " + Items[8]
      print string
      
      
