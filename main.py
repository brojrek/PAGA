from tkinter import ttk
from tkinter import filedialog as fd

import tkinter as tk
import datetime as dt
import time
import sys

from NW import *
import ncbiutils as ncbi
from Chromosome import *
from Population import *

#TODO Validate somewhere -> number of parents (got from selections - will be added in the future) 
# plus number of parents times a number of mutations must be LESS than population_size.
# i.e. if number of parents = 10
# and number of  mutations = 7
# population_size must be greater than 10 + 10*7 = 80 (in that case, 100 would be a good minimum population_size)
#TODO Input validation
#TODO Window displaying alignment progress live
#TODO Mutation split gap from section? Not sure

class Master:
    def __init__(self, master):
        self.master = master
        self.button_frame = ButtonsFrame(self)
        self.button_frame.button_file.configure(command = self.change_to_file)
        self.button_frame.button_search.configure(command = self.change_to_search)
        self.button_frame.button_search_id.configure(command = self.change_to_search_id)
        
        self.file_frame = FileFrame(self)
        self.search_frame = SearchFrame(self)
        self.search_id_frame = SearchIDFrame(self)

        self.active_sequence_frame = self.file_frame
        self.alignment_frame = AlignmentFrame(self)
        self.mutations_frame = MutationsFrame(self)

        self.frame = tk.Frame(self.master)

    def change_to_search(self):
        self.button_frame.button_search.configure(bg = "lightgreen")
        self.button_frame.button_search_id.configure(bg = "lightgrey")
        self.button_frame.button_file.configure(bg = "lightgrey")
        self.search_frame.frame.grid(row=0, column=0, padx=10, pady=90, rowspan=1)
        self.active_sequence_frame = self.search_frame

        self.file_frame.frame.grid_forget()
        self.search_id_frame.frame.grid_forget()
    
    def change_to_file(self):
        self.button_frame.button_file.configure(bg = "lightgreen")
        self.button_frame.button_search.configure(bg = "lightgrey")
        self.button_frame.button_search_id.configure(bg = "lightgrey")
        self.file_frame.frame.grid(row=0, column=0, padx=10, pady=90, rowspan=1)
        self.active_sequence_frame = self.file_frame

        self.search_frame.frame.grid_forget()
        self.search_id_frame.frame.grid_forget()

    def change_to_search_id(self):
        self.button_frame.button_file.configure(bg = "lightgrey")
        self.button_frame.button_search.configure(bg = "lightgrey")
        self.button_frame.button_search_id.configure(bg = "lightgreen")
        self.search_id_frame.frame.grid(row=0, column=0, padx=10, pady=90, rowspan=1)
        self.active_sequence_frame = self.search_id_frame

        self.file_frame.frame.grid_forget()
        self.search_frame.frame.grid_forget()
    
    def run_alignment(self):
        #TODO Validation of input parameters

        self.sequence_A = self.active_sequence_frame.sequence_A
        self.sequence_B = self.active_sequence_frame.sequence_B
        self.sequence_A_header = self.active_sequence_frame.sequence_A_header
        self.sequence_B_header = self.active_sequence_frame.sequence_B_header

        print("Active frame: " + str(type(self.active_sequence_frame)))
        print("Ran on sequences:")
        print(self.sequence_A)
        print(self.sequence_B)  
        if len(self.sequence_A) == 0:
            print("Something went wrong, no sequences passed to alignment.")
            return ""
        parent_chromosome = Chromosome(self.sequence_A, self.sequence_B, [self.alignment_frame.matrix, self.alignment_frame.alphabet], self.alignment_frame.gap_penalty.get())#, self.mutations_frame.gap_prolong_odds.get(), self.mutations_frame.gap_shuffle_odds.get(), self.mutations_frame.gap_remove_odds.get(), self.mutations_frame.gaps_limit_factor.get())
        print(parent_chromosome.score)

        # Set Chromosome proporties
        parent_chromosome.gap_prolong_odds = self.mutations_frame.gap_prolong_odds.get()
        parent_chromosome.new_gaps_limit_factor = self.mutations_frame.gaps_limit_factor.get()
        parent_chromosome.gap_shuffle_odds = self.mutations_frame.gap_shuffle_odds.get()
        parent_chromosome.gap_remove_odds = self.mutations_frame.gap_remove_odds.get()

        # Set Population and its proporties
        population = Population(parent_chromosome, self.alignment_frame.population_size.get())

        ## active mutations
        population.prolongation_active = self.mutations_frame.mutation_prolongation.get()
        population.add_one_gap_active = self.mutations_frame.mutation_add_single_gap.get()
        population.add_multiple_gaps_active = self.mutations_frame.mutation_add_multiple_gaps.get()
        population.gap_shuffle_gaps_active = self.mutations_frame.mutation_shuffle_gaps.get()
        population.move_gap_active = self.mutations_frame.mutation_move_gap.get()
        population.move_section_active = self.mutations_frame.mutation_move_gaps_section.get()
        population.gap_remove_active = self.mutations_frame.mutation_remove_gap.get()

        ## parent number
        population.parent_number = self.alignment_frame.parent_number.get()

        # Set stuff for genetic algorythm
        stop = False
        generations_without_improvement = 0
        highest_score = 0
        
        start_time = time.time()
        for n in range(self.alignment_frame.maximum_generations.get()):
            print('Iteration: ' + str(n) + ', ' + str(time.time() - start_time) + ' seconds passed. Highscore: ' + str(population.population[0].score))
            sys.stdout.flush()
            time.sleep(1)
            if highest_score < population.population[0].score:
                highest_score = population.population[0].score
                generations_without_improvement = 0
            else:
                generations_without_improvement += 1

            if generations_without_improvement > self.alignment_frame.max_generations_without_improvement.get():
                break

            if self.alignment_frame.selection.get() == "Tournament":
                population.new_generation(population.selection_tournament())
            elif self.alignment_frame.selection.get() == "Roulette":
                population.new_generation(population.selection_roulette())
            else: #self.alignment_frame.selection.get() == "Truncation"
                population.new_generation(population.selection_truncation())

        end_time = time.time()
        time_diff = round(end_time - start_time, 2)

        print("PAGA finished in " + str(n) + " generations.")
        print("Total execution time: " + str(time_diff) + " seconds.")
        #nw_results -> [[seq1, seq2], score]
        nw_results = nw_alignment(self.sequence_A, self.sequence_B, [self.alignment_frame.matrix, self.alignment_frame.alphabet], self.alignment_frame.gap_penalty.get())
        AlignmentResultWindow(self, [self.sequence_A_header, self.sequence_B_header], [self.sequence_A, self.sequence_B], time_diff, population.population[0], nw_results)

        
        print(population.population[0].score, population.population[0].sequence_A, population.population[0].sequence_B)
        #print([x.score for x in population.population])

class SearchFrame:
    def __init__(self, master):
        self.master = master.master
        self.frame = tk.Frame(self.master, width = self.master.winfo_width(), height = self.master.winfo_height())
        self.frame.grid(row=0, column=0, padx=10, pady=90, rowspan=1)

        self.sequence_A_id = tk.StringVar(value = "")
        self.sequence_B_id = tk.StringVar(value = "")
        self.sequence_A_header = tk.StringVar(value = "")
        self.sequence_B_header = tk.StringVar(value = "")

        self.sequence_A_header_display = tk.StringVar(value = "No sequence selected")
        self.sequence_B_header_display = tk.StringVar(value = "No sequence selected")

        self.sequence_A = ""
        self.sequence_B = ""

        #database
        self.database = tk.StringVar(value = "protein")
        database_label = tk.Label(self.frame, text = 'Database', font=('Times',11, 'bold'))
        database_entry = ttk.Combobox(
            self.frame,
            state = "readonly",
            values = ["protein", "nucleotide"], 
            font = ('Times 11'),
            width=17,
            textvariable = self.database
        )

        database_label.grid(row=1, column=0,  padx=10, pady=10, sticky="NW")
        database_entry.grid(row=1,column=1, pady=10, sticky="NE")

        #email
        self.email = tk.StringVar(value = "")
        email_label = tk.Label(self.frame, text = 'Email*', font=('Times', 11, 'bold'), padx=10, pady=10)
        email_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=self.email)

        email_label.grid(row=2,column=0, sticky="NW")
        email_entry.grid(row=2,column=1, pady=10, sticky="NE")

        #term
        term_label = tk.Label(self.frame, text = 'Term*', font=('Times', 11, 'bold'), padx=10, pady=10)
        term_entry = tk.Entry(self.frame, font = ('Times',11,'normal'))

        term_label.grid(row=3,column=0, sticky="NW")
        term_entry.grid(row=3,column=1, pady=10, sticky="NE")

        #organism
        organism_label = tk.Label(self.frame, text = 'Organism', font=('Times', 11, 'bold'), padx=10, pady=10)
        organism_entry = tk.Entry(self.frame, font = ('Times',11,'normal'))

        organism_label.grid(row=4,column=0, sticky="NW")
        organism_entry.grid(row=4,column=1, pady=10, sticky="NE")
        
        #minlen
        minlen = tk.IntVar(value = 0)
        minlen_label = tk.Label(self.frame, text = 'Min. length', font=('Times', 11, 'bold'), padx=10, pady=10)
        minlen_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=minlen)

        minlen_label.grid(row=5,column=0, sticky="NW")
        minlen_entry.grid(row=5,column=1, pady=10, sticky="NE")

        #maxlen
        maxlen = tk.IntVar(value = 999999)
        maxlen_label = tk.Label(self.frame, text = 'Max. length', font=('Times', 11, 'bold'), padx=10, pady=10)
        maxlen_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=maxlen)

        maxlen_label.grid(row=6,column=0, sticky="NW")
        maxlen_entry.grid(row=6,column=1, pady=10, sticky="NE")

        #fromdate
        fromdate = tk.StringVar(value = "1950/12/31")
        fromdate_label = tk.Label(self.frame, text = 'From [YYYY/MM/DD]', font=('Times', 11, 'bold'), padx=10, pady=10)
        fromdate_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=fromdate)

        fromdate_label.grid(row=7,column=0, sticky="NW")
        fromdate_entry.grid(row=7,column=1, pady=10, sticky="NE")

        #todate
        todate = tk.StringVar(value = str(dt.date.today().strftime('%Y/%m/%d')))
        todate_label = tk.Label(self.frame, text = 'To [YYYY/MM/DD]', font=('Times', 11, 'bold'), padx=10, pady=10)
        todate_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=todate)

        todate_label.grid(row=8,column=0, sticky="NW")
        todate_entry.grid(row=8,column=1, pady=10, sticky="NE")

        #retmax
        retmax = tk.IntVar(value = 30)
        retmax_label = tk.Label(self.frame, text = 'Max. results', font=('Times', 11, 'bold'), padx=10, pady=10)
        retmax_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=retmax)

        retmax_label.grid(row=9,column=0, sticky="NW")
        retmax_entry.grid(row=9,column=1, pady=10, sticky="NE")

        ## Buttons and sequences id

        # Sequence 1
        self.sequence_A_id = tk.StringVar(value = "")
        seq1_label = tk.Label(self.frame, text = 'Sequence 1. ID', font=('Times', 11, 'bold'))
        seq1_entry = tk.Entry(self.frame, textvariable = self.sequence_A_id, font = ('Times',11,'normal'), state="readonly")
        search_button1 = ttk.Button(self.frame, text='Search...', command=lambda : self.search_button_press(1, seq1_entry, email_entry.get(), database_entry.get(), term_entry.get(), organism_entry.get(), minlen_entry.get(), maxlen_entry.get(), fromdate_entry.get(), todate_entry.get(), retmax_entry.get()))

        seq1_label.grid(row=11,column=0, padx=10, pady=10, sticky="NW")
        seq1_entry.grid(row=11,column=1, pady=10, sticky="NW")
        search_button1.grid(row=11,column=2, pady=10, sticky="NE")

        #Sequence 2
        self.sequence_B_id = tk.StringVar(value = "")
        seq2_label = tk.Label(self.frame, text = 'Sequence 2. ID', font=('Times', 11, 'bold'))
        seq2_entry = tk.Entry(self.frame, textvariable = self.sequence_B_id, font = ('Times',11,'normal'), state="readonly")
        search_button2 = ttk.Button(self.frame, text='Search...', command=lambda : self.search_button_press(2, seq2_entry, email_entry.get(), database_entry.get(), term_entry.get(), organism_entry.get(), minlen_entry.get(), maxlen_entry.get(), fromdate_entry.get(), todate_entry.get(), retmax_entry.get()))

        seq2_label.grid(row=13,column=0, padx=10, pady=10, sticky="NW")
        seq2_entry.grid(row=13,column=1, pady=10, sticky="NW")
        search_button2.grid(row=13,column=2, pady=10, sticky="NE")

        #Loaded sequences
        loaded_A_label = tk.Label(self.frame, textvariable = self.sequence_A_header_display, font=('Times', 10), padx=10, pady=10)
        loaded_B_label = tk.Label(self.frame, textvariable = self.sequence_B_header_display, font=('Times', 10), padx=10, pady=10)

        loaded_A_label.grid(row=10,column=0, columnspan=3, sticky="W")
        loaded_B_label.grid(row=12,column=0, columnspan=3, sticky="W")

        self.frame.grid_forget()

    def display_search(self, results, sequence_number, entry):
        SearchResultWindow(self, results, sequence_number, entry)

    def search_button_press(self, sequence_number, entry, email, database, term, organism, minlen, maxlen, fromdate, todate, retmax):
        results = ncbi.search_button_clicked(email, database, term, organism, minlen, maxlen, fromdate, todate, retmax)
        self.display_search(results, sequence_number, entry)

    def selected_result(self, selected_id, selected_name, sequence_number, entry_id):
        
        if len(selected_name) > 50:
            selected_name_display = selected_name[:47] + "..."
        else:
            selected_name_display = selected_name

        if sequence_number == 1:
            self.sequence_A_id.set(value = str(selected_id))
            self.sequence_A_header.set(value = str(selected_name))
            self.sequence_A_header_display.set(value = str(selected_name_display))
            self.sequence_A = ncbi.sequence_from_id(str(selected_id), self.database.get(), self.email.get())[1]
            print(self.sequence_A)
        else:
            self.sequence_B_id.set(value = str(selected_id))
            self.sequence_B_header.set(value = str(selected_name))
            self.sequence_B_header_display.set(value = str(selected_name_display))
            self.sequence_B = ncbi.sequence_from_id(str(selected_id), self.database.get(), self.email.get())[1]
            print(self.sequence_B)
        entry_id.delete(0, "end")
        entry_id.insert(0, str(selected_id))

        print(selected_id, selected_name)

class SearchIDFrame:
    def __init__(self, master):
        self.master = master.master
        self.frame = tk.Frame(self.master, width = self.master.winfo_width(), height = self.master.winfo_height())
        self.frame.grid(row=0, column=0, padx=10, pady=90, rowspan=1)

        self.sequence_A = ""
        self.sequence_B = ""

        self.sequence_A_header = tk.StringVar(value = "")
        self.sequence_B_header = tk.StringVar(value = "")

        self.sequence_A_header_display = tk.StringVar(value = "No sequence selected")
        self.sequence_B_header_display = tk.StringVar(value = "No sequence selected")

        #database
        self.database = tk.StringVar(value="protein")
        database_label = tk.Label(self.frame, text = 'Database', font=('Times',11, 'bold'))
        database_entry = ttk.Combobox(
            self.frame,
            state = "readonly",
            values = ["protein", "nucleotide"], 
            font = ('Times 11'),
            width=17,
            textvariable=self.database
        )

        database_label.grid(row=1, column=0,  padx=10, pady=10, sticky="NW")
        database_entry.grid(row=1,column=1, pady=10, sticky="NE")

        #email
        self.email = tk.StringVar(value = "")
        email_label = tk.Label(self.frame, text = 'Email*', font=('Times', 11, 'bold'), padx=10, pady=10)
        email_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=self.email)

        email_label.grid(row=2,column=0, sticky="NW")
        email_entry.grid(row=2,column=1, pady=10, sticky="NE")

        #Sequence 1.
        self.sequence_A_id = tk.StringVar(value = "")

        loaded_A_label = tk.Label(self.frame, textvariable = self.sequence_A_header_display, font=('Times', 10), padx=10, pady=10)
        sequence_A_id_label = tk.Label(self.frame, text = 'Sequence 1. ID*', font=('Times', 11, 'bold'), padx=10, pady=10)
        sequence_A_id_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=self.sequence_A_id)
        search_buttonA = ttk.Button(self.frame, text='Search...', command=lambda : self.search_button_press(self.database.get(), self.email.get(), self.sequence_A_id.get(), 0))
        
        loaded_A_label.grid(row=3,column=0, columnspan=3, sticky="W")
        sequence_A_id_label.grid(row=4,column=0, sticky="NW")
        sequence_A_id_entry.grid(row=4,column=1, pady=10, sticky="NE")
        search_buttonA.grid(row=4,column=2, pady=10, sticky="NE")

        #Sequence 2.
        self.sequence_B_id = tk.StringVar(value = "")
        loaded_B_label = tk.Label(self.frame, textvariable = self.sequence_B_header_display, font=('Times', 10), padx=10, pady=10)
        sequence_B_id_label = tk.Label(self.frame, text = 'Sequence 2. ID*', font=('Times', 11, 'bold'), padx=10, pady=10)
        sequence_B_id_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=self.sequence_B_id)
        search_buttonB = ttk.Button(self.frame, text='Search...', command=lambda : self.search_button_press(self.database.get(), self.email.get(), self.sequence_B_id.get(), 1))
        
        loaded_B_label.grid(row=5,column=0, columnspan=3, sticky="W")
        sequence_B_id_label.grid(row=6,column=0, sticky="NW")
        sequence_B_id_entry.grid(row=6,column=1, pady=10, sticky="NE")
        search_buttonB.grid(row=6,column=2, pady=10, sticky="NE")

        self.frame.grid_forget() #Initially file frame will be displayed.
        
    def search_button_press(self, database, email, id, sequence_number):
        sequence_header, sequence = ncbi.sequence_from_id(id, database, email)
        print(sequence_header)
        print(sequence)
        if sequence_number == 0:
            self.sequence_A_header = tk.StringVar(value = sequence_header)
            if len(sequence_header) > 50:
                self.sequence_A_header_display.set(value = sequence_header[:47] + "...")
            else:
                self.sequence_A_header_display.set(value = sequence_header)
            self.sequence_A = sequence
            
        else: #sequence_number == 1:
            self.sequence_B_header = tk.StringVar(value = sequence_header)
            if len(sequence_header) > 50:
                self.sequence_B_header_display.set(value = sequence_header[:47] + "...")
            else:
                self.sequence_B_header_display.set(value = sequence_header)
            self.sequence_B = sequence

class FileFrame:
    def __init__(self, master):
        self.master = master.master
        self.frame = tk.Frame(self.master, width = self.master.winfo_width(), height = self.master.winfo_height())
        self.frame.grid(row=0, column=0, padx=10, pady=70, rowspan=1)


        self.sequence_A_id = tk.StringVar(value = "")
        self.sequence_B_id = tk.StringVar(value = "")

        self.sequence_A = ""
        self.sequence_B = ""

        self.sequence_A_header = tk.StringVar(value = "")
        self.sequence_B_header = tk.StringVar(value = "")

        self.sequence_A_header_display = tk.StringVar(value = "No sequence selected")
        self.sequence_B_header_display = tk.StringVar(value = "No sequence selected")


        #file1
        file1_label = tk.Label(self.frame, text = 'File 1', font=('Times', 11, 'bold'))
        file1_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), state='readonly', textvariable = self.sequence_A_id )
        file1_button = ttk.Button(self.frame, text='Open...', command=lambda : self.read_fasta_file(file1_entry, 1))

        #file2
        file2_label = tk.Label(self.frame, text = 'File 2', font=('Times', 11, 'bold'))
        file2_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), state='readonly', textvariable = self.sequence_B_id)
        file2_button = ttk.Button(self.frame, text='Open...', command=lambda : self.read_fasta_file(file2_entry, 2))

        file1_label.grid(row=2,column=0, padx=10, pady=10, sticky="NW")
        file1_entry.grid(row=2,column=1, pady=10, sticky="NW")
        file1_button.grid(row=2,column=2, pady=10, sticky="NE")

        file2_label.grid(row=4,column=0, padx=10, pady=10, sticky="NW")
        file2_entry.grid(row=4,column=1, pady=10, sticky="NW")
        file2_button.grid(row=4,column=2, pady=10, sticky="NE")

        #Loaded sequences
        loaded_A_label = tk.Label(self.frame, textvariable = self.sequence_A_header_display, font=('Times', 10), padx=10, pady=10)
        loaded_B_label = tk.Label(self.frame, textvariable = self.sequence_B_header_display, font=('Times', 10), padx=10, pady=10)

        loaded_A_label.grid(row=1,column=0, columnspan=3, sticky="W")
        loaded_B_label.grid(row=3,column=0, columnspan=3, sticky="W")


    def read_fasta_file(self, entry, n):
        filetypes = ( ('FASTA files', ['*.txt', '*.fasta']), ('All files', '*.*') )
        f = fd.askopenfile(filetypes=filetypes, initialdir="#Specify the file path")
        if self.read_fasta(f.name) == "":
            return ""
        
        if n == 1:
            name, self.sequence_A = self.read_fasta(f.name)
            self.sequence_A_header.set(name)
            if len(name) > 50:
                name = name[:47] + "..."
            self.sequence_A_header_display.set(name)
        else:
            name, self.sequence_B = self.read_fasta(f.name)
            self.sequence_B_header.set(name)
            if len(name) > 50:
                name = name[:47] + "..."
            self.sequence_B_header_display.set(name)

        entry.configure(state='normal')
        entry.delete(0, "end")
        entry.insert(0, f.name.split("/")[-1])
        entry.configure(state='readonly')

    def read_fasta(self, filename):
        '''
        Reads a sequence from .fasta file containing exactly one sequence. First line is a header.

        filename: string representing a path or name of the .fasta file

        returns: string representing header, string representing sequence read
        '''
        with open(filename, "r") as seq_file:
            f=(seq_file.read()).split('\n')
            seq=''

            if f[0][0] != ">": # Validation
                tk.messagebox.showerror(title = 'Error occured', message = 'Selected file\'s format is not FASTA.')
                return ""
            
            for i in range(1, len(f)):
                if f[i].count(">") > 0: # Validation
                    tk.messagebox.showerror(title = 'Error occured', message = 'Selected file contains more than one sequence.')
                    return ""
                seq=seq+f[i]

        return f[0], seq
        
class AlignmentFrame:
    def __init__(self, master):
        self.main = master
        self.master = master.master
        self.frame = tk.Frame(self.master, width = self.master.winfo_width(), height = self.master.winfo_height())
        self.frame.grid(row=0, column=1, padx=10, pady=10, rowspan=2, sticky="N")

        frame_alignment_label = tk.Label(self.frame, text = 'ALIGNMENT CONFIGURATION', font=('Times', 11, 'bold'), padx=10, pady=10)
        frame_alignment_label.grid(row=0,column=0, columnspan=2)

        #Population_size
        self.population_size = tk.IntVar(value = 100)
        population_size_label = tk.Label(self.frame, text = 'Population size', font=('Times', 11, 'bold'), padx=10, pady=10)
        population_size_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable = self.population_size)

        population_size_label.grid(row=1,column=0, sticky="W")
        population_size_entry.grid(row=1,column=1, sticky="W")

        #Maximum_generations
        self.maximum_generations = tk.IntVar(value = 100)
        maximum_generations_label = tk.Label(self.frame, text = 'Max. generations', font=('Times', 11, 'bold'), padx=10, pady=10)
        maximum_generations_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=self.maximum_generations)

        maximum_generations_label.grid(row=2,column=0, sticky="W")
        maximum_generations_entry.grid(row=2,column=1, sticky="W")

        #Max_generations_without_improvement
        self.max_generations_without_improvement = tk.IntVar(value = 20)
        max_generations_without_improvement_label = tk.Label(self.frame, text = 'Max. w/o improvement', font=('Times', 11, 'bold'), padx=10, pady=10)
        max_generations_without_improvement_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable=self.max_generations_without_improvement)

        max_generations_without_improvement_label.grid(row=3,column=0, sticky="W")
        max_generations_without_improvement_entry.grid(row=3,column=1, sticky="W")

        #Selection
        self.selection = tk.StringVar(value = "Tournament")
        selection_label = tk.Label(self.frame, text = 'Selection method', font=('Times',11, 'bold'), padx=10, pady=10)
        selection_entry = ttk.Combobox(
            self.frame,
            state = "readonly",
            textvariable = self.selection,
            values = ["Tournament", "Truncation", "Roulette"], 
            font = ('Times 11'),
            width=18
        )
        #selection_entry.set("Tournament")

        selection_label.grid(row=4,column=0, sticky="W")
        selection_entry.grid(row=4,column=1, sticky="W")

        #Parent_number
        self.parent_number = tk.IntVar(value = 5)
        parents_number_label = tk.Label(self.frame, text = 'Parents selected', font=('Times', 11, 'bold'), padx=10, pady=10)
        parents_number_entry = tk.Entry(self.frame, font = ('Times',11,'normal'), textvariable = self.parent_number)

        parents_number_label.grid(row=5,column=0, sticky="W")
        parents_number_entry.grid(row=5,column=1, sticky="W")

        #Gap_penalty
        self.gap_penalty = tk.IntVar(value = 5)
        gap_penalty_label = tk.Label(self.frame, text = 'Gap penalty', font=('Times', 11, 'bold'), padx=10, pady=10)
        gap_penalty_entry = tk.Entry(self.frame, textvariable = self.gap_penalty, font = ('Times',11,'normal'))

        gap_penalty_label.grid(row=6,column=0, sticky="W")
        gap_penalty_entry.grid(row=6,column=1, sticky="W")

        #Matrix
        self.matrix_filename = tk.StringVar(value = "")
        matrix_filename_label = tk.Label(self.frame, text = 'Score matrix', font=('Times', 11, 'bold'), padx=10, pady=10)
        matrix_filename_entry = tk.Entry(self.frame, textvariable = self.matrix_filename, font = ('Times',11,'normal'), state = 'readonly')

        matrix_filename_label.grid(row=7,column=0, sticky="W")
        matrix_filename_entry.grid(row=7,column=1, sticky="W")

        matrix_button = ttk.Button(self.frame, text='Open...', command = lambda: self.read_matrix_file(matrix_filename_entry))
        matrix_button.grid(row=7,column=2, pady=10, sticky="NE")

        #Run
        run_style = ttk.Style()
        run_style.configure('my.TButton', font=('Times', 20))
        run_button = ttk.Button(self.frame, text='Align!', width=20, style='my.TButton', command = master.run_alignment)
        run_button.grid(row=10,column=0, pady=10,columnspan=2, sticky="NE")

    def read_matrix_file(self, entry):
        '''
        #TODO documentation
        '''
        filetypes = ( ('Text files', '*.txt'), ('All files', '*.*') )
        f = fd.askopenfile(filetypes=filetypes, initialdir="#Specify the file path")
        if self.read_matrix(f.name) == "":
            return ""

        self.matrix, self.alphabet = self.read_matrix(f.name)

        entry.configure(state='normal')
        entry.delete(0, "end")
        entry.insert(0, f.name.split("/")[-1])
        entry.configure(state='readonly')

    def read_matrix(self, filename):
        '''
        #TODO documentation
        '''
        with open(filename, "r") as f:
            content = f.read()
            content = content.splitlines()
            n = len(content)-1
            matrix=np.zeros((n,n),'i') #d - float i - int  
            f.close()
            alphabet={}
            row=content[0].split()
            for j in range(0,len(row)):
                alphabet[row[j]]=j

            for i in range(1,len(content)):
                row=content[i].split()
                for j in range(1,len(row)):
                    if row[j]!='':
                        matrix[i-1,j-1] = float(row[j])
                        matrix[j-1,i-1] = float(matrix[i-1,j-1])

        return matrix, alphabet


class MutationsFrame:
    def __init__(self, master):
        self.master = master.master

        self.frame = tk.Frame(self.master, width = self.master.winfo_width(), height = self.master.winfo_height())
        self.frame.grid(row=0, column=2, padx=10, pady=10, rowspan=2, sticky="N")

        self.frame_label = tk.Label(self.frame, text = 'ACTIVE MUTATIONS', font=('Times', 11, 'bold'), padx=10, pady=10)
        self.frame_label.grid(row=0, column=0)

        # Prolongation
        self.mutation_prolongation = tk.BooleanVar(value=True)
        self.c_prolongation = tk.Checkbutton(self.frame, text='Prolong existing gaps',variable=self.mutation_prolongation, onvalue=True, offvalue=False)
        self.c_prolongation.grid(row = 1, column = 0, sticky="W")
        
        #gap_prolong_odds
        self.gap_prolong_odds = tk.IntVar(value = 6)
        gap_prolong_odds_label = tk.Label(self.frame, text = 'Gap prolong odds', font=('Times', 10, 'normal'), padx=10, pady=10)
        gap_prolong_odds_entry = tk.Entry(self.frame, font = ('Times',10,'normal'), width = 5, textvariable=self.gap_prolong_odds)

        gap_prolong_odds_label.grid(row=2,column=0, sticky="W")
        gap_prolong_odds_entry.grid(row=2,column=0, sticky="E")
          
        # Add one gap
        self.mutation_add_single_gap = tk.BooleanVar(value=True)
        c_add_single_gap = tk.Checkbutton(self.frame, text='Add single gap',variable=self.mutation_add_single_gap, onvalue=True, offvalue=False)
        c_add_single_gap.grid(row = 3, column = 0, sticky="W")

        # Add multiple gaps
        self.mutation_add_multiple_gaps = tk.BooleanVar(value=True)
        c_add_multiple_gap = tk.Checkbutton(self.frame, text='Add multiple gaps',variable=self.mutation_add_multiple_gaps, onvalue=True, offvalue=False)
        c_add_multiple_gap.grid(row = 4, column = 0, sticky="W")
        
        #new_gaps_limit_factor
        self.gaps_limit_factor = tk.IntVar(value = 10)
        new_gaps_limit_factor_label = tk.Label(self.frame, text = 'New gaps limit factor', font=('Times', 10, 'normal'), padx=10, pady=10)
        new_gaps_limit_factor_entry = tk.Entry(self.frame, font = ('Times',10,'normal'), width = 5, textvariable=self.gaps_limit_factor)

        new_gaps_limit_factor_label.grid(row=5,column=0, sticky="W")
        new_gaps_limit_factor_entry.grid(row=5,column=0, sticky="E")

        # Shuffle gaps
        self.mutation_shuffle_gaps = tk.BooleanVar(value=True)
        c_shuffle_gaps = tk.Checkbutton(self.frame, text='Shuffle existing gaps',variable=self.mutation_shuffle_gaps, onvalue=True, offvalue=False)
        c_shuffle_gaps.grid(row = 6, column = 0, sticky="W")
        
        #gap_shuffle_odds
        self.gap_shuffle_odds = tk.IntVar(value = 3)
        gap_shuffle_odds_label = tk.Label(self.frame, text = 'Gap shuffle odds', font=('Times', 10, 'normal'), padx=10, pady=10)
        gap_shuffle_odds_entry = tk.Entry(self.frame, font = ('Times',10,'normal'), width = 5, textvariable = self.gap_shuffle_odds)

        gap_shuffle_odds_label.grid(row=7,column=0, sticky="W")
        gap_shuffle_odds_entry.grid(row=7,column=0, sticky="E")
        
        # Move gap
        self.mutation_move_gap = tk.BooleanVar(value=True)
        c_move_gap = tk.Checkbutton(self.frame, text='Move single gap',variable=self.mutation_move_gap, onvalue=True, offvalue=False)
        c_move_gap.grid(row = 8, column = 0, sticky="W")

        # Move a section of gaps
        self.mutation_move_gaps_section = tk.BooleanVar(value=True)
        c_move_gaps_section = tk.Checkbutton(self.frame, text='Move a section of gaps',variable = self.mutation_move_gaps_section, onvalue=True, offvalue=False)
        c_move_gaps_section.grid(row = 9, column = 0, sticky="W")

        # Remove gaps
        self.mutation_remove_gap = tk.BooleanVar(value=True)
        c_remove_gap = tk.Checkbutton(self.frame, text='Remove gaps',variable=self.mutation_remove_gap, onvalue=True, offvalue=False)
        c_remove_gap.grid(row = 10, column = 0, sticky="W")
        
        #gap_remove_odds
        self.gap_remove_odds = tk.IntVar(value = 4)
        gap_remove_odds_label = tk.Label(self.frame, text = 'Gap remove odds', font=('Times', 10, 'normal'), padx=10, pady=10)
        gap_remove_odds_entry = tk.Entry(self.frame, font = ('Times',10,'normal'), width = 5, textvariable=self.gap_remove_odds)

        gap_remove_odds_label.grid(row=11,column=0, sticky="W")
        gap_remove_odds_entry.grid(row=11,column=0, sticky="E")

class ButtonsFrame:
    def __init__(self, master):
        self.master = master.master

        self.frame = tk.Frame(self.master, width = self.master.winfo_width(), height = self.master.winfo_height())
        self.frame.grid(row=0, column=0, padx=10, pady=20, sticky="N")

        buttons_label = tk.Label(self.frame, text="SOURCE OF SEQUENCES FOR ALIGNMENT", font=('Times', 11, 'bold'), padx=10, pady=0)

        #FILE
        self.button_file = tk.Button(self.frame, text = 'FILE', font=('Times', 11, 'bold'), padx=1, pady=0, bg ="lightgreen")

        #SEARCH
        self.button_search = tk.Button(self.frame, text = 'SEARCH', font=('Times', 11, 'bold'), padx=1, pady=0)
        
        #SEARCH_ID
        self.button_search_id = tk.Button(self.frame, text = 'SEARCH ID', font=('Times', 11, 'bold'), padx=1, pady=0)

        buttons_label.grid(row=0, column=0, columnspan = 2)
        self.button_file.grid(row=1,column=0, sticky="w")
        self.button_search.grid(row=1,column=0, sticky="E")
        self.button_search_id.grid(row=1,column=1, sticky="W")

class SearchResultWindow:
    def __init__(self, master, results, sequence_number, entry):
        self.search = master
        self.master = master.master
        self.sequence_number = sequence_number
        self.entry = entry

        self.result_window = tk.Toplevel(self.master)
 
        self.result_window.title('Search results')
 
        # sets the geometry of toplevel
        self.result_window.geometry('600x350')

        # define columns
        columns = (' ', 'ID', 'Name', 'Length')
        self.tree = ttk.Treeview(self.result_window, columns=columns, show='headings', selectmode='browse')
        self.tree['columns'] = columns

        # format columns
        self.tree.column("#0", width=0,  stretch='no')
        self.tree.column(" ",anchor='center', width=15)
        self.tree.column("ID",anchor='center', width=110)
        self.tree.column("Name",anchor='center',width=360)
        self.tree.column("Length",anchor='center',width=100)
        #self.tree.column("Data dodania",anchor='center',width=150)
        
        # headers
        self.tree.heading("#0",text="",anchor='center')
        self.tree.heading(" ", text = "", anchor='center')
        self.tree.heading("ID",text="ID",anchor='center')
        self.tree.heading("Name",text="Name",anchor='center')
        self.tree.heading("Length",text="Length",anchor='center')
 
        #add data
        for n in range(0, len(results)):
            #print(len(results[n]))
            self.tree.insert(parent = '', index = 'end', iid = n, text = '', values = (n+1, results[n][0], results[n][1], len(results[n][2])))

        
        scrollbar = ttk.Scrollbar(self.result_window, command=self.tree.yview)
        self.tree.configure(yscroll=scrollbar.set)
        scrollbar.grid(column=4, row=0, sticky='ns')
        self.tree.bind('<Double-Button-1>', self.select_item)

        self.tree.grid(column=0, row=0, sticky='news')

    def select_item(self, event):
        current_item = self.tree.item(self.tree.focus())
        print(current_item)
        self.search.selected_result(current_item['values'][1], current_item['values'][2], self.sequence_number, self.entry)
        
        self.result_window.destroy()
        self.result_window.update()

class AlignmentResultWindow:
    def __init__(self, master, input_sequences_headers, input_sequences, ga_time, ga_results_chromosome, nw_results):
        self.master = master.master

        self.result_window = tk.Toplevel(self.master)
 
        self.result_window.title('Alignment results')
    
        self.result_window.geometry("800x410")
        self.result_window.resizable('True', 'True')

        self.frame = tk.Frame(self.result_window)
        self.frame.grid()
        
        #scrollbar
        sb = tk.Scrollbar(self.frame, orient = 'horizontal')
        #sb.pack(side = 'bottom', fill = 'x')
        sb.grid(row = 1, column = 0, sticky = tk.EW)
        
        # nw_results -> [[seq1, seq2], score]
        text0 = str("Input headers:\nSeq1: " + input_sequences_headers[0].get() + "\nSeq2: " + input_sequences_headers[1].get() + "\n\n")
        text1 = str("Input sequences:\nSeq1: " + input_sequences[0] + "\nSeq2: " + input_sequences[1])
        text2 = str('\n\nNW results:\nScore: ' + str(nw_results[1]) + '\nSeq1: ' + str(nw_results[0][0]) + '\nSeq2: ' + str(nw_results[0][1]) + '\n\n')
        text3 = str('Genetic algorithm results:\nScore: ' + str(ga_results_chromosome.score) + '\nSeq1: ' + str(ga_results_chromosome.sequence_A) + '\nSeq2: ' + str(ga_results_chromosome.sequence_B))
        text4 = str('\n\nGenetic algorithm execution time: ' + str(ga_time) + ' seconds.')
        text5 = str("\nUsed score matrix (file name): " + str(master.alignment_frame.matrix_filename.get()))
        self.textALL = text0 + text1 + text2 + text3 + text4 + text5
        
        text = tk.Text(self.frame, font="TkFixedFont", wrap = 'none', xscrollcommand=sb.set)
        
        text.grid(row=0, column=0)
        sb.config(command = text.xview)

        text.insert('end', self.textALL)
        
        #Save button
        self.button_save = ttk.Button(self.frame, text='Save as...', command = lambda: self.button_save_pressed(self.textALL))

        self.button_save.grid(row=0,column=1, sticky="S")

        tk.mainloop()

    def button_save_pressed(self, text):
        # ask filename to be used
        filename = fd.asksaveasfilename(defaultextension="txt", initialfile="alignment_" + str(dt.date.today().strftime("%m%d")))
        if filename:
            with open(filename, mode = 'w', encoding = 'UTF-8') as file:
                file.write(text)


def main():
    root = tk.Tk()
    root.title("PAGA - Pairwise Alignment by Genetic Algortihm")
    window_width = 1030
    window_height = 640

    # get the screen dimension
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    if screen_width > screen_width:
        window_width = screen_width

    if screen_height > screen_height:
        window_height = screen_height

    # find the center point
    center_x = int(screen_width/2 - window_width/2)
    center_y = int(screen_height/2 - window_height/2)

    # set the position of the window to the center of the screen
    root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')
    root.resizable('False', 'False')
    root.grid_columnconfigure(0, weight=1)
    root.grid_columnconfigure(1, weight=1)
    root.grid_columnconfigure(2, weight=1)

    app = Master(root)
    root.mainloop()

if __name__ == '__main__':
    main()

