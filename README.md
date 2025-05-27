# ICA
introduce group project ICA -- Genetic Sequence Analysis and Mutation Simulation Tool   poster   coding
This program requires us to analysis any given RNA sequence that starts with AUG.
Firstly, we need to identify the most frequent trinucleotide and report the most frequent DNA trinucleotide.
Secondly, find the encoded amino acid for the most frequent trinucleotide based on the first function.
Then, draw the bar chart of the amino acid frequency
And finally, as for the additional function, we simulate random point mutations to assess their effects on RNA translation.

This code first imports the Counter class from Python's collections module, which is specifically designed for counting. We create a list that contains some numbers. Then, we use Counter to convert this list into a dictionary, where the keys are the elements in the list and the values are the counts of these elements.
Next, we print this dictionary, which tells us the count of each number. Then, we call the most_common() method, which returns a list of elements sorted in descending order by their counts. If we pass a parameter to most_common(), such as most_common(1), it will only return the element with the highest count.
Finally, we access this list by index to obtain the element with the highest count and its count.
Such tools can help us better understand and process information in data analysis.