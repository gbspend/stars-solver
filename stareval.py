import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from copy import deepcopy
from heapq import heappush, heappop

X = "X"
O = "%"

#===============================================================================

def draw(grid, labels, title=None):
    rows = len(grid)
    cols = len(grid[0])

    if len(labels) != rows or any(len(row) != cols for row in labels):
        raise ValueError("grid and labels must have the same dimensions")
    
    print()
    if title: print(title)

    for r in range(rows):

        # --- Draw horizontal boundaries above row r ---
        line = "+"
        for c in range(cols):
            if r == 0 or grid[r][c] != grid[r-1][c]:
                line += "---+"
            else:
                line += "   +"
        print(line)

        # --- Draw cell contents with vertical boundaries ---
        line = ""
        for c in range(cols):
            if c == 0 or grid[r][c] != grid[r][c-1]:
                line += "|"
            else:
                line += " "

            label = labels[r][c] if labels[r][c] else " "
            line += f" {label} "

        line += "|"  # right border
        print(line)

    # --- Bottom border ---
    line = "+"
    for c in range(cols):
        line += "---+"
    print(line)

#===============================================================================

def neighbors(row, col, max_size):
    neighbors = []
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            # Skip the cell itself
            if dr == 0 and dc == 0:
                continue
            new_r = row + dr
            new_c = col + dc
            # Check bounds
            if 0 <= new_r < max_size and 0 <= new_c < max_size:
                neighbors.append((new_r, new_c))
    return neighbors

def empty_labels(grid):
    ret = []
    for _ in grid:
        curr = []
        for __ in _:
            curr.append("")
        ret.append(curr)
    return ret

def neigh_set(grid):
    arr_g=[]
    for a in grid:
        arr_g+=a
    return set(arr_g)

#NOT in Puzzle class because we'll be trying with temp labels
#check grid for any errors (not whether it's complete)
def check_errors(grid, labels):
    MAX = len(grid)
    areas = defaultdict(int)
    cols = defaultdict(int)
    #first loop checks touching stars and >2 per r/c/a
    for r in range(MAX):
        row_count = 0
        for c in range(MAX):
            if labels[r][c] == O:
                row_count += 1
                cols[c] += 1
                areas[grid[r][c]] += 1
                if row_count > 2 or cols[c] > 2 or areas[c] > 2:
                    return False
                for nr, nc in neighbors(r,c,MAX):
                    if labels[nr][nc] == O:
                        return False
                        
    #second loop checks if blanks + stars < 2
    areas = defaultdict(list) #list all empty spots
    area_stars = defaultdict(int)
    cols = defaultdict(list)
    col_stars = defaultdict(int)
    for r in range(MAX):
        row_spots = []
        row_stars = 0
        for c in range(MAX):
            if labels[r][c] == O:
                row_stars += 1
                col_stars[c] += 1
                area_stars[grid[r][c]] += 1
            if labels[r][c] == "":
                t = (r,c)
                row_spots.append(t)
                cols[c].append(t)
                areas[grid[r][c]].append(t)
            if r == MAX-1: #on last row, we're finished with the col
                if not sufficient_clumps(grid, labels,cols[c],col_stars[c]):
                    return False
        if not sufficient_clumps(grid, labels,row_spots,row_stars):
            return False
    for n in neigh_set(grid): #CHECK ALL AREAS, ESP. IF NOT THERE!
        if not sufficient_clumps(grid, labels,areas[n],area_stars[n]):
            return False
    return True

#assumes check_errors is True (i.e. no errors)
#marks ONE or TWO stars in same r/c/a (if possible) then returns True
#   if no stars found, returns False
#if a r/c (maybe area later?) has 2 stars, it fills in the rest x
#run it while return >0
def find_next_star(grid, labels):
    MAX = len(grid)
    areas = defaultdict(list) #list all empty spots
    area_stars = defaultdict(int)
    cols = defaultdict(list)
    col_stars = defaultdict(int)
    for r in range(MAX):
        row_spots = []
        row_stars = 0
        for c in range(MAX):
            if labels[r][c] == O:
                row_stars += 1
                col_stars[c] += 1
                area_stars[grid[r][c]] += 1
            if labels[r][c] == "":
                t = (r,c)
                row_spots.append(t)
                cols[c].append(t)
                areas[grid[r][c]].append(t)
            if r == MAX-1: #on last row, we're finished with the col
                num = check_line_clumps(labels,cols[c],col_stars[c])
                if num > 0:
                    return num
        num = check_line_clumps(labels,row_spots,row_stars)
        if num > 0:
            return num
    for n in areas:
        num = check_area_clumps(grid, labels,areas[n],area_stars[n])
        if num > 0:
            return num
    return 0

#helper for check_line_clumps
#autodetect r/c by which changes \+1
def contig(curr, old):
    if curr[0] == old[0] and curr[1] == old[1] + 1:
        return True
    if curr[1] == old[1] and curr[0] == old[0] + 1:
        return True
    return False

#helpers for both "check" and "sufficient" funcs
def get_line_clumps(labels, spots):
    old = None
    clumps = [] #list of clumps
    clump = [] #current clump
    for i in range(len(spots)):
        curr = spots[i]
        if old is None or contig(curr, old):
            clump.append(curr)
        else:
            clumps.append(clump)
            clump = [curr]
        old = curr
    clumps.append(clump)
    return clumps

def get_area_clumps(labels, spots):
    clumps = [] #list of clumps
    for i in range(len(spots)):
        curr = spots[i]
        foundClump = False
        r,c = curr
        for n in neighbors(r,c,len(labels)):
            for cl in clumps:
                if n in cl:
                    cl.append(curr)
                    foundClump = True
                    break
            if foundClump:
                break
        if not foundClump:
            clumps.append([curr]) #add "seed clump"
    return clumps

#these are easy...ish...
#modifies lablels in place
def check_line_clumps(labels, spots, stars):
    if stars > 2:
        print("CHECK_LINE_CLUMPS found stars > 3!")
        exit(1)
    if stars == 2: #fill in xs
        for r,c in spots:
            labels[r][c] = X
        return 0
    #so stars is 0-1 here
    missing_stars = 2-stars
    ret = 0
    clumps = get_line_clumps(labels, spots)
    if len(clumps) == missing_stars:
        allGood = True
        for clump in clumps:
            if len(clump) > 2:
                allGood = False
        if allGood:
            for clump in clumps:
                if len(clump) == 1:
                    add_star(labels, clump[0])
                    ret = 1
    #is this overly specific (and I'm missing something)
    #   or actually a very specific smallest case?
    if missing_stars == 2 and len(clumps) == 1 and len(clumps[0]) == 3:
        add_star(labels, clumps[0][0])
        add_star(labels, clumps[0][-1])
        ret = 2
    return ret

#...hmmm area clumps are tricky...
#somethingsomething 3x3 area...?
def check_area_clumps(grid,labels, spots, stars):
    if stars > 2:
        print("CHECK_AREA_CLUMPS found stars > 3!") 
        exit(1)
    if stars == 2: #fill in xs
        for r,c in spots:
            labels[r][c] = X
        return 0
    #so stars is 0-1 here
    missing_stars = 2-stars
    ret = 0
    clumps = get_area_clumps(labels, spots)
    if len(clumps) == missing_stars:
        allGood = True
        for clump in clumps:
            if len(clump) > 2:
                allGood = False
        if allGood:
            for clump in clumps:
                if len(clump) == 1:
                    add_star(labels, clump[0])
                    ret+=1
    return ret

def sufficient_clumps(grid, labels, spots, stars):
    if stars > 2:
        return False
    if stars == 2:
        return True
    if stars == 1:
        return len(spots) > 0
    else: #stars == 0
        for curr in spots:
            r,c = curr
            neighs = neighbors(r,c,len(labels))
            rem = set(neighs + [curr])
            left = set(spots) - rem
            if len(left) > 0:
                return True
        return False
    
#returns true if no neighboring stars, false otherwise (avoid)
def add_star(labels, t):
    r,c = t
    neighs = neighbors(r,c,len(labels))
    #check first
    for nr, nc in neighs:
        if labels[nr][nc] == O:
            return False
    for nr, nc in neighs:
        labels[nr][nc] = X
    labels[r][c] = O
    return True

#assumes check_errors is True (i.e. no errors)
#finds as many as possible, but still run it while returns >0
#returns number of xs found
def find_xs(grid, labels):
    MAX = len(grid)
    xs = 0
    for r in range(MAX):
        for c in range(MAX):
            if labels[r][c] == "":
                cp = deepcopy(labels)
                add_star(cp, (r,c))
                curr = (r,c)
                if not check_errors(grid, cp):
                    labels[r][c] = X
                    xs += 1
    return xs

#requires no errors
#returns True if no errors, else returns # of stars
#   use "if check_complete(...) is True", not ==
def check_complete(grid, labels):
    MAX = len(grid)
    total = 0
    for r in range(MAX):
        for c in range(MAX):
            if labels[r][c] == O:
                total += 1
    if total == MAX * 2:
        return True
    else:
        return total
    

class StarEval:
    def __init__(self):
        self.reset()
    
    def reset(self):
        self.star_seq = []
        self.x_seq = []
        self.uncert_count = 0
        self.max_depth = 0
        self.child_hist = [] #record different uncert/recurs results...?
        self.complete = False
    
    def print_report(self):
        print(self.complete)
        print(self.star_seq)
        print(self.x_seq)
        print(self.uncert_count, self.max_depth)
        if self.child_hist:
            print("Children history:")
            for ch in self.child_hist:
                print(ch)
    
    def update_depth(self, i):
        if i > self.max_depth:
            self.max_depth = i
            
    def note_children(self, ch):
        self.child_hist.append([-t[0] for t in ch])
        
    def uncert_search(self, grid, labels, curr_stars, i=0):
        self.uncert_count+=1
        self.update_depth(i)
        MAX = len(grid)
        children = []
        for r in range(MAX):
            for c in range(MAX):
                if labels[r][c] == "":
                    cp = deepcopy(labels)
                    add_star(cp, (r,c))
                    #cp should NEVER have errors here
                    cp_labels, num = self.solve(grid, cp, False)
                    if num is True: #SOLVED!
                        self.note_children(children)
                        return cp_labels
                    elif num > curr_stars+1: #+1 because we forced one!
                        heappush(children, (-num, cp_labels)) 
        #can peek with children[0], maybe return something more than None...
        if children:
            self.note_children(children)
            child_stars, child_labels = heappop(children)
            child_stars *= -1 #undo max heap
            result = self.uncert_search(grid, child_labels, child_stars, i+1)
            if result:
                return result
        return None
        
    def solve(self, grid, labels=None, deep=True):
        if deep: self.reset()
        if labels is None:
            labels = empty_labels(grid)
        while True:
            changed = False
            ss = 0 #find stars can only do 1 or 2; so aggregate
            while True:
                temp = find_next_star(grid, labels)
                if temp > 0:
                    ss += temp
                    if not check_errors(grid,labels):
                        if deep:
                            print("ERRORS IN SOLVE LOOP (after stars)!")
                            draw(grid,labels)
                            exit(1) #this is a bug IF this is the first time (unforced)
                        else: #we forced a star, and it caused an error: OK, keep searching!
                            return None, 0
                    changed = True
                else:
                    if deep:
                        self.star_seq.append(ss)
                    break
            while True:
                xs = find_xs(grid,labels)
                if xs > 0:
                    if not check_errors(grid,labels):
                        if deep: 
                            print("ERRORS IN SOLVE LOOP (after Xs)!")
                            draw(grid,labels)
                            exit(1) #this is a bug IF this is the first time (unforced)
                        else: #we forced a star, and it caused an error: OK, keep searching!
                            return None, 0
                    elif deep:
                        self.x_seq.append(xs)
                    changed = True
                else:
                    break
            if not changed:
                break
        comp = check_complete(grid, labels) #comp is True if done, # stars otherwise
        if comp is True:
            if deep: draw(grid,labels, "FINISHED")
            self.complete = True
        else:
            if deep: draw(grid,labels, "INCOMPLETE")
            if deep:
                labels = self.uncert_search(grid,labels,comp)
                if labels:
                    if not check_errors(grid, labels):
                        draw(grid,labels,"ERRORS IN SEARCH-PSEUDO-SOLVED!!")
                        exit(1)
                    if not check_complete(grid, labels):
                        draw(grid,labels,"INCOMPLETE SEARCH-PSEUDO-SOLVED!!")
                        exit(1)
                    draw(grid,labels, "SEARCH-SOLVED!!")
                    return labels, True
                else:
                    print("DEEPER SEARCH COULDN'T SOLVE :(")
                    return None, False
            else:
                pass#draw(grid,labels, "INCOMPLETE")
        return labels, comp

#NOTES
#Binary space partition to split non, hand-drawn regions!
#https://en.wikipedia.org/wiki/Binary_space_partitioning
