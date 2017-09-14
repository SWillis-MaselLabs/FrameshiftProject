# The purpose of this submodule is to take in a list that contains multiple sublists, each of length
# 2, and to combine the lists based on whether they have a trivial intersection or not:
# Example:
# input list = [[a,b],[c,d],[c,e],[e,f]]
# output list = [[a,b], [c,d,e,f]]
# Each entry in each sublist should only appear once in both the sublist of the output
# and only once in the sublist it is contained in

def SetConsolidator(l):
    # SetConsolidator takes in a list that has multiple sublists, each of length 2
    # It defines output to be an empty list that it will use to merge lists that share a common
    # element, and add new lists that have a trivial intersection with the rest of the sublists
    output = []
    # Because output starts out blank, we want to be able to use a variable to tell it that
    # nothing has been added yet. 
    firstItem = True

    # we start out by examining the items in the input list "l"
    for item in l:
        # we declare a variable that will be able to keep track of whether an element
        # of any sublist has been appended to the output list. This way we avoid adding an element or
        # sublist multiple times. The numberAppended variable is updated to 0 at the beginning of
        # each loop through the sublists of l
        Appended = False

        # Because we want to be able to compare sublists, we sort each sublist so it's in a
        # standardized format
        sortedItem = sorted(item)

        # The elements of each sublist are the UIDs for two genes. We declare variables
        # so we don't have to use indexed sublists throughout the code.
        gene1 = sortedItem[0]
        gene2 = sortedItem[1]

        # If nothing has been appended to the output list yet, then we add the first sublist
        # of l to the output and set the firstItem variable to False so that nothing more will
        # be added without examination first. 
        if firstItem == True:
            output.append(sortedItem)
            firstItem = False
        # If it's not the first thing added to list, then we need to look at what is already in the
        # output before we add anything more to it.
        else:
            # Before we start looking through each individual sublist in the output, we can do a
            # quick check to see if the item coming from l is already in the output. If it is, we
            # discard it and get a new item from l.
            if sortedItem in output:
                pass
            # if the sorted item is not in the output list, then we look through all the
            
            else:
            # To examine the output list, we loop through the output list to see what's in it
                for sublist in output:
                    # If gene1 is in the sublist already, we see if gene2 is in it as well
                    if gene1 in sublist:
                        # if gene 2 isn't in the sublist, then we add it and make a note that
                        # the sublist had something appended to the overall list
                        if gene2 not in sublist:
                            sublist.append(gene2)
                            Appended = True
                        # if gene2 is also in the sublist, we make a note that both genes
                        # were found in the output list
                        else:
                            Appended = True
                    # if gene1 is not in the sublist, we check to see if gene2 is in the sublist
                    else:
                        # if it is, we add gene 1
                        if gene2 in sublist:
                            sublist.append(gene1)
                            Appended = True

                # the last thing we do before we look at the next sublist is see whether anything was added
                # to the output. If nothing was added to the output, then we add this sublist to the
                # output. We do this until we have exhausted all options
                if Appended == False:
                    output.append(sortedItem)

    # we can't return the output yet, though, because there's some funny business that may have
    # come about. Consider the following situation:

    # first sublist = [A,B]
    # second sublist = [C,D]
    # third sublist = [C,B]

    # with the above algorithm, if we looked at these lists in this order, we'd wind up
    # with multiple lists with a nontrivial intersection: [A,B,C] and [C,D,B] We need to
    # go through the output list and check to see if any such occurance took place

    # to do this, we check each sublist in the output list
    for sublist in output:
        # and using a second for loop, we cross-reference that sublist with all other
        # sublists in the output
        for sublist2 in output:
            # the first sublist will always have a non-trivial intersection with itself,
            # so we want to discard these cases
            if sublist == sublist2:
                pass
            # if the list is not comparing itself to itself, then we check to see the size
            # of the intersection
            else:
                intersection = list(set(sublist) & set(sublist2))
                print(len(intersection))
                # if the intersection is trivial, we move on
                if len(intersection) == 0:
                    pass
                # otherwise, we need to consolidate
                else:
                    # To consolidate, we consider every element in the second sublist and
                    # see if it's in the first sublist
                    for element in sublist2:
                        # if the element in sublist2 is not in the first sublist, then we add
                        # it to the first sublist
                        if element not in sublist:
                            sublist.append(element)
                        # if it is in the first sublist (and therefore in the intersection), then
                        # we move on
                        else:
                            pass
                    # once every element that wasn't in the first sublist that was in the second sublist
                    # is added, we remove sublist2 from the output. By doing this repeatedly, we remove
                    # all sublists that contain duplicate entries and wind up with the desired
                    # output
                    output.remove(sublist2)
    # finally we return the desired result
    return output
                        
