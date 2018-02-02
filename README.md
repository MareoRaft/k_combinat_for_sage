The goal is to translate Morse's code from (Maple?) to Sage.

Please read /Applications/sage/sage-7.2/src/sage/interfaces/maple.py.
You can also type `maple.fsolve('x^2=cos(x)+4', 'x=0..5')` into the Sage shell to get an explanation of how you should install/config Maple on your system.  One question is whether this maple integration works only in the console, or in any Sage script.


It may be possible to import maple into Sage, which would solve the problem at least for personal use!



Plan of attack
-------------------
  1. obtain code
  2. obtain Maple
  3. create Maple program that runs each function against certain inputs (of your choice, i guess) and stores the output to a friendly format (JSON).
  4. Use the JSON to generate a python pytest file
  5. translate the code, using the test file to make sure it works

Translation
-----------------
There are many possibilities.  One could translate by hand, which would allow for code optimization / programmer improvements.

One could build up a list of useful regexs in the process to translate common patterns to their Sage equivalents.

One could build a Maple parser to translate to Sage, but I'm worried this would be too hard :(

Whatever we do, we should document the code as we go.

