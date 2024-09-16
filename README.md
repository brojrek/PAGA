## Table of contents
* [About](#about)
* [Technologies](#technologies)
* [How to use](#how-to-use)

## General info
PAGA (Pairwise Alignment with Genetic Algorithm) is an implementation of genetic algorithm for sequence alignment.

The purpose of this implementation is of explorational character. While there exist tools performing pairwise alignment relatively fast and returning optimal result (i.e. Needleman-Wunsh algorithm), there is not such "perfect" tool for multiple sequences alignment.

According to https://pubmed.ncbi.nlm.nih.gov/11747615/, multiple sequence alignment is a NP-hard problem, therefore finding a heuristic solution (like a genetic algorithm) to this problem is a reasonable approach. PAGA's purpose is to propose and verify this approach for a pair of sequences with possibility to further develop a multiple sequences alignment use.

## Technologies
PAGA is written in Python 3.12.2 and uses threading, biopython and numpy libraries. GUI is created using tkinter.

## How to use

Launch main.py to start.

![Alt text](/README_images/Main_window.png?raw=true "Main window image")

### Buttons
Frame 1 contains three buttons. Pressing any of them will change the source of sequences for alignment depending on it's description.

### Search frame
Frame 2 is responsible for sequences input for alignment. Displayed frame depends on button selected in frame 1.

#### File source
Default sequences source is file. Frame 2 displays two fields allowing to load a sequences for alignment from selected file.
File must be of FASTA format containing exactly one sequence. Upon succesfully loading a file, header will be displayed above the corresponding field and file name will be put in that field. Sequence will be stored in memory in of frame 2 (file source).

#### Search source
Upon pressing "SEARCH" button in frame 1, frame 2 will display fields responsible for building a search query for NCBI database and
two fields allowing to search (via "Search.." button) and load selected results to frame 2 (search source) memory. Loaded sequence header will be displayed above corresponding field.

Fields with an asterisk symbol (*) are mandatory since required by NCBI search (Bio.Entrez). It is recommended to provide a valid email address.

#### Search ID source
Upon pressing "SEARCH ID" button in frame 1, frame 2 will display fields responsible for a search input in NCBI database and
two fields allowing to search by ID (via "Search.." button) and load selected results to frame 2 (search id source) memory. Loaded sequence header will be displayed above corresponding field.

Fields with an asterisk symbol (*) are mandatory since required by NCBI search (Bio.Entrez). It is recommended to provide a valid email address.

### Alignment frame
Frame 3 is responsible for alignment settings:
    - population size : a total number of chromosomes population
    - max generations : a maximum number of generations in genetic algorithm. If stop condition won't be triggered, this is the number of generations after which algorithm will stop and return results
    - max w/o improvement : a maximum number of consecutive generations without improvement of highscore after which algorithm will stop and return results
    - selection method : selection method used by genetic algorithm. Tournament is recommended choice.
    - parents selected : number of chromosomes selected by selection method to initiate next generation. Every generation is is initiated with selected parents. Then every parent is subjected to every active mutation and result (offspring) is appended to population
    - gap penalty : a positive number, that is subtracted from chromosome's score for every gap present at any index.
    - score matrix : loaded from a file with "Open..." button. For proper work it must contain every character present in chosen sequences (that is, all aminoacids for protein sequences, all nucleotides for nucleotide sequence) and scores for every pair

### Mutations frame
Frame 4 is reponsible for selecting and setting active mutations:
    - prolong existing gaps : insert a new gap near another, already present gap with given odds. Each gap in a sequence will be prolonged with odds of 1:(gap prolong odds)
    - add single gap : insert a new gap ar random index in a randomly selected sequnece in a chromosome
    - add multiple gaps : insert a random number of gaps on random indexes in a randomly selected sequence in a chromosome. Value of "new gaps limit factor" defines the upper limit of number of gaps that will be added as follows: ((length of a sequence) / (new gaps limit factor)) + 1
    - shuffle existing gaps : shuffle randomly selected gaps in randomly selected sequence in a chromosome to random indexes.
    Each gap in a sequence will be shuffled with odds of 1:(gap shuffle odds)
    - move single gap : randomly selects a single gap and moves it by one index in random direction
    - move a section of gaps : randomly selects a section of gaps and moves it by one index in random direction
    - remove gaps : removes a random number of gaps from a sequence. Each gap in a sequence will be removed with odds of 1:(gap remove odds)

### Align! button
Align! button runs an alignment. As input sequences it uses these loaded in currently active source frame. Both sequences and alignment settings, including scoring matrix, must be provided.

After alignment is finished, result window will show up. You can save results to .txt file using "Save as..." button.

### Tips
- You can also set a really (really, like 9999999) maximum number of generations and reasonable number (like 100) of maximum generations without improvement. Algorithm will eventually stop returning local maximum it found, not limited by maximum number of generations. This will be especially helpful with particularly demanding alignments of long and diverse sequences.

For example, aligning sequences with ids 2462630655 and 2272752337 with maximum generations = 1000 and maximum generations without improvement = 40 resulted with highscore of 1000 out of 8525 (according to NW alignment). Same alignment with maximum generations = 300000 and maximum generations without improvement = 100 resulted in score of 8390 out of 8525 (in 347 iterations, 28 seconds).

Don't be afraid to set ridiculously high values of these parameters!

- One of the flaws of the genetic algorithm is that it may focus too much on a single local maximum, which may be far from global maximum. Try to run the same alignment multiple times looking for better results.

- Progress bar showed during alignment represents a current generation number divided by maximum number of generations. Usually it won't fill to 100% since alignment will stop after set number of generations without improvement. Also, if you set a very high maximum generations number, you may not notice it progressing at all.
