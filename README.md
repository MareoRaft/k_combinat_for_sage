The goal is to translate Morse's code from (Maple?) to Sage.

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

