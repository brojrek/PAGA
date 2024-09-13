from Bio import Entrez
from tkinter import messagebox
from tkinter import ttk
from datetime import datetime
import tkinter as tk
import sys
from threading import Thread

def sequence_from_id(ID, database, email):
    Entrez.email = email
    searchResultHandle = Entrez.esearch(db = database, term = ID, retmax = 1)  #retmax=1000)
    searchResult = Entrez.read(searchResultHandle)
    try:
        handle = Entrez.efetch(db = database, id = searchResult["IdList"][0], rettype="fasta", retmode="text")
    except IndexError:
        err = 'Sequence with given ID not found, ID: ' + ID
        messagebox.showinfo('Message', err)
        return ''

    if ID == searchResult["IdList"][0]:
        print("Found a sequence")
        
    else:
        err = 'Sequence with given ID not found, ID: ' + ID
        messagebox.showinfo('Message', err)
        return ''
    
    seq = handle.read().split("\n", 1)
    return seq[0], seq[-1].replace('\n', '')

class SearchWindow:
    def __init__(self, master):
        self.master = master
        self.results = []
        self.progressbar = ttk.Progressbar(master = master, orient=tk.HORIZONTAL, length=160, mode="indeterminate")
        self.progressbar.grid(row = 1, padx = 5)

        label = ttk.Label(master, text="Searching...", font=20)
        label.grid(row = 0)
        self.progressbar.start()

    def search(self, email, database, term, organism, minlen, maxlen, fromdate, todate, retmax):

        # Search for sequences in a new thread
        t = Thread(target= lambda: self.searchseq(email, database, term, organism, minlen, maxlen, fromdate, todate, retmax))
        t.start()

        # Start checking periodically if the thread has finished.
        self.schedule_check(t)
    
    def schedule_check(self, t):
        """
        Schedule the execution of the `check_if_done()` function after
        one second.
        """
        self.master.after(1000, lambda: self.check_if_done(t))

    def check_if_done(self, t):
        # If the thread has finished, close the window
        if not t.is_alive():
            self.progressbar.stop()
            self.master.withdraw()
            self.master.quit()
            self.master.update() 
        else:
            # Otherwise check again after one second.
            self.schedule_check(t)

    def searchseq(self, email, database, term, organism, minlen, maxlen, fromdate, todate, retmax):
        email = str(email)
        if len(email) < 3 or '@' not in email:
            messagebox.showinfo("Message", "Enter a correct email.")
            return None
        database = str(database)
        if database != 'protein' and database != 'nucleotide':
            messagebox.showinfo("Message", "Select a correct database.")
            return None
        searchterm = str(term)
        if len(searchterm) < 1:
            messagebox.showinfo("Message", "Term is a mandatory field and must be filled.")
            return None
        organism = str(organism)
        # if:

        
        #filtr = str(entry_filtr.get())
        # if:

        
        #keyword = str(entry_keyword.get())
        # if:

        
        minlen = str(minlen)
        if minlen.isnumeric() == False:
            messagebox.showinfo("Message", "Given value is incorrect. Minimal length set to 0.")
            minlen = '0'
        maxlen = str(maxlen)
        if maxlen.isnumeric() == False:
            messagebox.showinfo("Message", "Given value is incorrect. Maximum length set to 99999999.")
            maxlen = '99999999'
        fromdate = str(fromdate)
        if len(fromdate) != 10 and len(fromdate) != 0: #TODO better valdiation
            messagebox.showinfo("Message", "Given value is incorrect. Start date will be ommited.")
            fromdate = ''

        todate = str(todate)
        
        if len(str(todate)) != 10 and todate.replace('/','').isnumeric() == False:
            messagebox.showinfo("Message", "Given value is incorrect. End date will be ommited.")
            todate = str(datetime.today().strftime('%Y/%m/%d'))
        
        retmax = str(retmax)
        if retmax.isnumeric() == False:
            messagebox.showinfo("Message", "Given value is incorrect. Max results set to 20.")
            retmax = '20'

        # Validation finished
        
        Entrez.email = email
        
        if len(organism) > 0:
            searchterm = searchterm + ' AND ' + organism + '[Organism]'

        #if len(filtr) > 0:
        #    searchterm = searchterm + ' AND ' + filtr + '[Filter]'

        #if len(keyword) > 0:
        #    searchterm = searchterm + ' AND ' + keyword + '[Keyword]'
            
        # min and max length.
        if minlen != '0' or maxlen != '123456789':
            searchterm = searchterm + ' AND ' + minlen + '[Sequence Length]:' + maxlen + '[Sequence Length]'
        
        # startdate and todate (YYYY/MM/DD)
        if fromdate != '' or todate != str(datetime.today().strftime('%Y/%m/%d')):
            searchterm = searchterm + ' AND ' + fromdate + '[Publication Date]:' + todate + '[Publication Date]'

        #SEARCHDB
        print(searchterm)
        searchResultHandle = Entrez.esearch(db = database, term = searchterm, retmax = int(retmax))  #retmax=1000)
        searchResult = Entrez.read(searchResultHandle)
        ids = searchResult["IdList"]
        results = []
        for id in ids:
            handle = Entrez.efetch(db = database, id=id, rettype="fasta", retmode="text")
            seq = handle.read().split("\n", 1)
            results.append([id, seq[0], seq[1].replace('\n', '')])
            handle.close()
        self.results = results

def search_button_clicked(email, database, term, organism, minlen, maxlen, fromdate, todate, retmax):
    root = tk.Tk()
    root.title("NCBI database search")
    window_width = 170
    window_height = 60

    # get the screen dimension
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # find the center point
    center_x = int(screen_width/2 - window_width/2)
    center_y = int(screen_height/2 - window_height/2)

    # set the position of the window to the center of the screen
    root.geometry(f'170x60+{center_x}+{center_y}')
    root.resizable('False', 'False')

    app = SearchWindow(root)
    app.search(email, database, term, organism, minlen, maxlen, fromdate, todate, retmax)

    root.mainloop()

    results = app.results

    root.destroy()
    return results
