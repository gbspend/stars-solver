import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from copy import deepcopy
from heapq import heappush, heappop

X = "×"
O = "☆"

class Node:
    def __init__(self, grid, labels):
        self.grid = grid
        self.labels = labels
        self.children = [] #expand very conservatively!
        
# ==================================================================

def draw(grid, labels, title=""):
    """
    grid   : 2D array of integers representing contiguous regions
    labels : 2D array of strings (same size as grid), 
             each element is '' or a single-character string
    """
    grid = np.array(grid)
    labels = np.array(labels)

    rows, cols = grid.shape

    if labels.shape != grid.shape:
        raise ValueError("grid and labels must have the same dimensions")

    fc = None
    if title == "FINISHED":
        fc = "#ffff99"
    fig, ax = plt.subplots(facecolor=fc)

    # Draw thin black grid lines
    for r in range(rows + 1):
        ax.plot([0, cols], [r, r], linewidth=0.5, color='black')
    for c in range(cols + 1):
        ax.plot([c, c], [0, rows], linewidth=0.5, color='black')

    # Draw thick black region boundaries
    for r in range(rows):
        for c in range(cols):
            current = grid[r, c]

            # Right boundary
            if c < cols - 1 and grid[r, c + 1] != current:
                ax.plot([c + 1, c + 1], [r, r + 1], linewidth=3, color='black')

            # Bottom boundary
            if r < rows - 1 and grid[r + 1, c] != current:
                ax.plot([c, c + 1], [r + 1, r + 1], linewidth=3, color='black')

    # Outer border (thick black)
    ax.plot([0, cols], [0, 0], linewidth=3, color='black')
    ax.plot([0, cols], [rows, rows], linewidth=3, color='black')
    ax.plot([0, 0], [0, rows], linewidth=3, color='black')
    ax.plot([cols, cols], [0, rows], linewidth=3, color='black')

    # Draw labels
    for r in range(rows):
        for c in range(cols):
            text = labels[r, c]
            if text:  # nonempty string
                ax.text(
                    c + 0.5,
                    r + 0.5,
                    text,
                    ha='center',
                    va='center',
                    fontsize=14,
                    color='black'
                )

    ax.set_xlim(0, cols)
    ax.set_ylim(rows, 0)
    ax.set_aspect('equal')
    ax.axis('off')
    
    if title:
        plt.title(title)
    plt.show()

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
                if not sufficient_clumps(labels,cols[c],col_stars[c]):
                    return False
        if not sufficient_clumps(labels,row_spots,row_stars):
            return False
    for n in neigh_set(grid): #CHECK ALL AREAS, ESP. IF NOT THERE!
        if not sufficient_clumps(labels,areas[n],area_stars[n]):
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
                if check_line_clumps(labels,cols[c],col_stars[c]):
                    return True
        if check_line_clumps(labels,row_spots,row_stars):
            return True
    for n in areas:
        if check_area_clumps(grid, labels,areas[n],area_stars[n]):
            return True
    return False

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
        return False
    #so stars is 0-1 here
    missing_stars = 2-stars
    ret = False
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
                    ret = True
    #is this overly specific (and I'm missing something)
    #   or actually a very specific smallest case?
    if missing_stars == 2 and len(clumps) == 1 and len(clumps[0]) == 3:
        add_star(labels, clumps[0][0])
        add_star(labels, clumps[0][-1])
        ret = True
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
        return False
    #so stars is 0-1 here
    missing_stars = 2-stars
    ret = False
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
                    ret = True

def sufficient_clumps(labels, spots, stars):
    if stars > 2:
        print("SUFF_CLUMPS found stars > 3!") 
        exit(1)
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

def solve1(grid):
    labels = empty_labels(grid)
    while True:
        changed = False
        while find_xs(grid,labels) > 0:
            if not check_errors(grid,labels):
                print("ERRORS IN SOLVE LOOP (after Xs)!")
                draw(grid,labels)
                exit(1)
            changed = True
        while find_next_star(grid, labels) > 0:
            if not check_errors(grid,labels):
                print("ERRORS IN SOLVE LOOP (after stars)!")
                draw(grid,labels)
                exit(1)
            changed = True
        if not changed:
            break
    draw(grid,labels, "FINISHED")
    return labels

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
    
def uncert_search(grid, labels, curr_stars, i=0):
    MAX = len(grid)
    children = []
    for r in range(MAX):
        for c in range(MAX):
            if labels[r][c] == "":
                cp = deepcopy(labels)
                add_star(cp, (r,c))
                #cp should NEVER have errors here
                cp_labels, num = solve2(grid, cp, False)
                if num is True: #SOLVED!
                    return cp_labels
                elif num > curr_stars+1: #+1 because we forced one!
                    heappush(children, (num, cp_labels))
                #use Node.children as heapq! Push on with num as key, ONLY if num > curr_stars
    #can peek with children[0], maybe return something more than None...
    if children:
        child_stars, child_labels = heappop(children)
        result = uncert_search(child_labels, child_stars, i+1)
        if result:
            return result
    return None

def solve2(grid, labels=None, deep=True):
    if labels is None:
        labels = empty_labels(grid)
    while True:
        changed = False
        while find_xs(grid,labels) > 0:
            if not deep:
                pass#draw(grid,labels, "FOUND X IN SEARCH")
            if not check_errors(grid,labels):
                if deep: 
                    print("ERRORS IN SOLVE LOOP (after Xs)!")
                    draw(grid,labels)
                    exit(1) #this is a bug IF this is the first time (unforced)
                else: #we forced a star, and it caused an error: OK, keep searching!
                    return None, 0
            changed = True
        while find_next_star(grid, labels) > 0:
            if not deep:
                pass#draw(grid,labels, "FOUND STAR IN SEARCH")
            if not check_errors(grid,labels):
                if deep:
                    print("ERRORS IN SOLVE LOOP (after stars)!")
                    draw(grid,labels)
                    exit(1) #this is a bug IF this is the first time (unforced)
                else: #we forced a star, and it caused an error: OK, keep searching!
                    return None, 0
            changed = True
        if not changed:
            break
    comp = check_complete(grid, labels) #comp is True if done, # stars otherwise
    if comp is True:
        if deep: draw(grid,labels, "FINISHED")
    else:
        if deep: draw(grid,labels, "INCOMPLETE")
        if deep:
            labels = uncert_search(grid,labels,comp)
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

#===============================================================================

if __name__ == "__main__":
    test1 = [
        [1,1,2,2,2,2,3,3,3],
        [1,1,1,1,1,1,3,3,3],
        [4,4,4,4,1,1,1,1,1],
        [1,1,1,1,1,5,5,6,1],
        [1,7,7,7,7,5,5,6,1],
        [1,1,1,1,1,5,5,6,1],
        [1,8,8,1,1,1,1,6,1],
        [1,8,8,1,9,9,9,9,1],
        [1,8,8,1,1,1,1,1,1]
    ]
    test1_2 = [
        [1,1,1,1,2,2,2,2,2],
        [1,3,3,1,1,2,4,4,2],
        [1,1,3,3,1,2,4,4,2],
        [1,5,1,1,1,2,4,4,2],
        [5,5,5,5,5,5,2,2,2],
        [6,6,6,6,5,5,7,7,7],
        [6,6,8,6,6,7,7,9,9],
        [6,8,8,8,6,7,9,9,7],
        [6,8,8,8,6,7,7,7,7],
    ]
    test2 = [
        [1,1,1,1,1,5,5,5,5,6],
        [1,2,2,2,1,7,7,7,5,6],
        [1,2,2,2,1,7,7,6,5,6],
        [1,2,2,2,1,6,6,6,5,6],
        [1,1,1,1,1,6,6,6,6,6],
        [8,8,8,9,9,3,3,3,3,3],
        [10,10,8,8,9,3,4,4,4,3],
        [10,10,9,8,9,3,4,4,4,3],
        [10,9,9,8,9,3,4,4,4,3],
        [9,9,9,9,9,3,3,3,3,3],
    ]
    test3 = [
        [1,1,1,1,2,3,3,3,3],
        [1,4,4,4,2,2,2,2,3],
        [1,4,4,4,2,2,2,2,3],
        [1,4,4,4,5,5,2,2,3],
        [4,4,4,5,5,6,6,6,6],
        [7,8,8,8,8,6,6,6,9],
        [7,8,8,8,8,6,6,6,9],
        [7,8,8,8,8,6,6,6,9],
        [7,7,7,7,8,9,9,9,9],
    ]
    test4 = [
        [1,1,1,1,1,1,1,1,1,1],
        [1,4,4,4,5,5,6,6,6,1],
        [1,7,4,4,5,5,6,6,6,1],
        [1,7,7,4,5,8,8,6,6,1],
        [1,7,7,4,5,8,8,6,6,1],
        [1,7,3,8,5,8,8,0,0,1],
        [1,7,3,8,8,8,8,0,9,9],
        [1,3,3,3,3,2,0,0,9,9],
        [1,3,2,2,2,2,0,0,9,9],
        [2,2,2,2,2,2,2,9,9,9],
    ]
    #l1 = empty_labels(test1)
    solve2(test4)
    #l2 = empty_labels(test2)
    #draw(test2, l2)
    #l2 = solve1(test2)
    #draw(test2, l2)
    print("bon weekend!")


#TODO: add branching tests!
#   remember: sort of BFS... like pursue branch with FEWEST possibilities first