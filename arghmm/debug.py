



def add_arg_thread2(arg, new_name, thread, recombs, arg3=None):


    def is_local_coal(arg, node, pos, local):
        return (len(node.children) == 2 and
                node.children[0] in local and
                arg.get_local_parent(node.children[0], pos-.5) == node and
                node.children[1] in local and
                arg.get_local_parent(node.children[1], pos-.5) == node and
                node.children[0] != node.children[1])



    def walk_up(arg, leaves, time, pos, ignore=None):

        print
        print "walk_up", leaves, time, ignore

        order = dict((node, i) for i, node in enumerate(
            arg.postorder_marginal_tree(pos-.5)))
        local = set(order.keys())
        if ignore is not None and ignore in arg:
            ptr = arg[ignore]
        else:
            ptr = None
        if ptr in local:
            local.remove(ptr)
            ptr = arg.get_local_parent(ptr, pos-.5)
            
            while ptr and ptr in local:
                if (len(ptr.children) == 2 and
                    ((ptr.children[0] in local and
                      arg.get_local_parent(ptr.children[0], pos-.5) == ptr) or
                     (ptr.children[1] in local and
                      arg.get_local_parent(ptr.children[1], pos-.5) == ptr))):
                    print " halt", ptr
                    break
                print " remove", ptr
                local.remove(ptr)
                ptr = arg.get_local_parent(ptr, pos-.5)

        print " local", local

        queue = [(order[arg[x]], arg[x]) for x in leaves]
        seen = set(x[1] for x in queue)
        heapq.heapify(queue)

        while len(queue) > 1:
            print "queue", queue
            i, node = heapq.heappop(queue)
            parent = arg.get_local_parent(node, pos-.5)
            if parent and parent not in seen:
                seen.add(parent)
                heapq.heappush(queue, (order[parent], parent))
        node = queue[0][1]
        parent = arg.get_local_parent(node, pos-.5)

        print " node", node, node.age
        assert node.age <= time
        
        while parent and parent.age <= time:
            if is_local_coal(arg, parent, pos, local):
                print " stop", node.age, parent.age, parent.children
                break            
            node = parent
            parent = arg.get_local_parent(node, pos-.5)
            print " node", node, node.age

        if parent:
            if parent.age < time:
                print leaves, parent.age, time, ignore
                tree = arg.get_marginal_tree(pos-.5).get_tree()
                tree.write()
                treelib.draw_tree_names(tree, maxlen=8, minlen=8)
                assert False

        return node


    def add_node(arg, node, time, pos, event):

        assert node.age <= time, (node.age, time)

        node2 = arg.new_node(event=event, age=time, children=[node], pos=pos)
        if event == "coal":
            node2.pos = 0

        parent = arg.get_local_parent(node, pos-.5)
        if parent:
            assert time <= parent.age, (time, parent.age)
            node.parents[node.parents.index(parent)] = node2
            parent.children[parent.children.index(node)] = node2
            node2.parents.append(parent)
        else:
            node.parents.append(node2)

        return node2


    arg_recomb = dict((x.pos, x) for x in iter_visible_recombs(arg))
    recomb_clades = [
        (pos-1, None) + get_clade_point(arg, rnode, rtime, pos-1)
        for pos, rnode, rtime in recombs] + [
        (node.pos, node.name) +
        get_clade_point(arg, node.name, node.age, node.pos)
        for node in iter_visible_recombs(arg)]
    recomb_clades.sort()

    # make initial tree
    arg2 = arg.get_marginal_tree(-1)
    arglib.remove_single_lineages(arg2)

    start = get_clade_point(arg, thread[0][0], thread[0][1], 0)
    node = walk_up(arg2, start[0], start[1], -1)
    node2 = add_node(arg2, node, start[1], -1, "coal")
    leaf = arg2.new_node(name=new_name, event="gene", age=0)
    leaf.parents.append(node2)
    node2.children.append(leaf)

    print "init arg2"
    tree = arg2.get_marginal_tree(-.5).get_tree()
    phylo.hash_order_tree(tree)
    tree.write()

    if arg3:
        print "init arg3"
        tree = arg3.get_marginal_tree(-.5).get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        tree.write()


    # ensure all thread changes occur at recomb points
    r = set(x[0] for x in recomb_clades)
    for pos in range(1, len(thread)):
        if thread[pos] != thread[pos-1]:
            assert (pos-1) in r, (pos, thread[pos], thread[pos-1], sorted(r))


    # add each recomb and re-coal
    for rpos, rname, rleaves, rtime in recomb_clades:
        print "------------------------------------------"
        print "recomb=", (rpos, rleaves, rtime), (rpos in arg_recomb)
        print "thread=", thread[rpos], thread[rpos+1]

        for node in arg2:
            if node.event == "recomb":
                assert len(node.parents) == 2, node
        
        if rpos in arg_recomb:
            # find re-coal for existing recomb

            if thread[rpos][1] != thread[rpos+1][1]:
                if rtime > min(thread[rpos][1], thread[rpos+1][1]):
                    print ">>", rtime, thread[rpos], thread[rpos+1]
                    treelib.draw_tree_names(
                        arg.get_marginal_tree(rpos-.5).get_tree(),
                        maxlen=8, minlen=8)
                    treelib.draw_tree_names(
                        arg.get_marginal_tree(rpos+.5).get_tree(),
                    maxlen=8, minlen=8)
                    assert False
            
            node = arg_recomb[rpos]
            local2 = set(arg.postorder_marginal_tree(rpos+.5))
            last = node
            node = arg.get_local_parent(node, rpos+.5)
            while (not is_local_coal(arg, node, rpos+1, local2)):
                last = node
                node = arg.get_local_parent(node, rpos+.5)
            c = node.children
            child = c[0] if c[1] == last else c[1]
            recoal = node

            print ">>", node, c
            treelib.draw_tree_names(
                arg.get_marginal_tree(rpos-.5).get_tree(),
                maxlen=8, minlen=8)
            treelib.draw_tree_names(
                arg.get_marginal_tree(rpos+.5).get_tree(),
                maxlen=8, minlen=8)
                
            cleaves, ctime = get_clade_point(
                arg, child.name, node.age, rpos-.5)

            # ensure this is not a cycle
            assert set(rleaves) != set(cleaves), rleaves

            # get local tree T^{n-1}_i and add new branch
            tree = arg.get_marginal_tree(rpos+.5)
            arglib.remove_single_lineages(tree)            
            node_name, time = thread[rpos+1]
            node = tree[node_name]

            node2 = add_node(tree, node, time, rpos+1, "coal")
            if not node2.parents:
                tree.root = node2
            leaf = tree.new_node(name=new_name, event="gene", age=0)
            leaf.parents.append(node2)
            node2.children.append(leaf)

            print "tmp", (rleaves, rtime), (cleaves, ctime), thread[rpos+1]
            tree2 = tree.get_tree()
            phylo.hash_order_tree(tree2)
            tree2.write()

            recomb = walk_up(tree, rleaves, rtime, rpos+1, new_name)

            if recomb == node2 and rtime == node2.age:
                # recomb and new coal-state are near each other
                # we must decide if recomb goes above or below coal-state

                # if this is a mediated SPR, then recomb goes below.
                # otherwise it goes above.

                # SPR is mediated if previous coal state is not recomb branch
                
                print "recoal=", recoal.name, "node2=", node2.name
                treelib.draw_tree_names(
                    arg.get_marginal_tree(rpos+.5).get_tree(),
                    maxlen=8, minlen=8)
                treelib.draw_tree_names(tree2, maxlen=8, minlen=8)

                node_name, time = thread[rpos]

                print "thread", thread[rpos], thread[rpos+1]
                
                if node2.children[0].name != node_name:
                    # this is a mediated coal
                    recomb = node2.children[0]

            
            coal = recomb.parents[0]
            c = coal.children
            child = c[0] if c[1] == recomb else c[1]


            print ">", (list(tree.leaf_names(recomb)), rtime), \
                  (list(tree.leaf_names(child)), coal.age)
            
            # get coal point in T^n_i
            rleaves, rtime = get_clade_point(
                tree, recomb.name, rtime, rpos+1)
            cleaves, ctime = get_clade_point(
                tree, child.name, coal.age, rpos+1)

            print ">>> arg2"
            tree = arg2.get_marginal_tree(rpos+.5).get_tree()
            phylo.hash_order_tree(tree)
            treelib.draw_tree_names(tree, minlen=8, maxlen=8)
            treelib.remove_single_children(tree)
            tree.write()

            node1 = walk_up(arg2, rleaves, rtime, rpos+1)
            node2 = walk_up(arg2, cleaves, ctime, rpos+1, node1.name)

    
        else:
            # find re-coal for new recomb
            
            assert rtime <= thread[rpos][1], (rtime, thread[rpos][1])
            
            if rleaves == [new_name]:
                # recomb on new branch, coal given thread
                cleaves, ctime = get_clade_point(
                    arg, thread[rpos+1][0], thread[rpos+1][1], rpos+.5)
                assert ctime >= rtime, (rtime, ctime)
                
                node1 = walk_up(arg2, rleaves, rtime, rpos+1)
                node2 = walk_up(arg2, cleaves, ctime, rpos+1, new_name)
                
            else:
                # recomb in ARG, coal on new branch
                cleaves = [new_name]
                ctime = thread[rpos+1][1]
                assert ctime >= rtime, (rtime, ctime)

                # NOTE: new_name is not ignored for walk_up on rleaves
                # because I do not want the recombination to be higher
                # than the coal point, which could happen if the recomb time
                # is the same as the current coal time.
                node1 = walk_up(arg2, rleaves, rtime, rpos+1)
                node2 = walk_up(arg2, cleaves, ctime, rpos+1, node1.name)


        print "add", rpos, rpos in arg_recomb, (rleaves, rtime), (cleaves, ctime)
        print "  node1", list(arg2.leaf_names(node1))
        print "  node2", list(arg2.leaf_names(node2))


        if arg3:
            print "arg3"
            tree = arg3.get_marginal_tree(rpos+.5).get_tree()
            treelib.remove_single_children(tree)
            phylo.hash_order_tree(tree)
            tree.write()

        print "arg"
        tree = arg.get_marginal_tree(rpos+.5).get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        tree.write()

        print "arg (i-1)"
        tree = arg.get_marginal_tree(rpos-.5).get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        tree.write()

        #if set(rleaves) == set(cleaves):
        #    print "skip", rleaves
        #    continue

        assert node1.parents
        assert rtime <= ctime

        recomb = add_node(arg2, node1, rtime, rpos, "recomb")
        if node1 == node2:
            node2 = recomb
        coal = add_node(arg2, node2, ctime, rpos, "coal")

        recomb.parents.append(coal)
        coal.children.append(recomb)

        #arglib.assert_arg(arg2)

        print "arg2"
        tree = arg2.get_marginal_tree(rpos+.5).get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        tree.write()

        print "arg2 (i-1)"
        tree = arg2.get_marginal_tree(rpos-.5).get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        tree.write()


        print "  r", recomb, recomb.children, recomb.parents
        print "  c", coal, coal.children, coal.parents


        node, time = get_coal_point(arg2, arg2[new_name], rpos+1)
        assert time == thread[rpos+1][1], (time, thread[rpos+1][1])

    
    
    return arg2
    


def calc_transition_probs2(tree, states, nlineages, times,
                           time_steps, popsizes, rho):

    ntimes = len(time_steps)
    #treelen = sum(x.get_dist() for x in tree)
    treelen = get_treelen(tree, times)
    mintime = time_steps[0]
    nbranches, nrecombs, ncoals = nlineages
    
    # A_{k,j} =& s'_{j-2} k_{j-2} / (2N) + \sum_{m=k}^{j-3} s'_m k_m / (2N) 
    #         =& s'_{j-2} k_{j-2} / (2N) + A_{k,j-1}.
    
    A = util.make_matrix(ntimes, ntimes, 0.0)
    for k in xrange(ntimes):
        # A[k][k] = A[k][k+1] = 0
        for j in xrange(k+2, ntimes):
            l = j - 2
            A[k][j] = A[k][j-1] + time_steps[l] * nbranches[l] / (2.0 * popsizes[l])

    # B_{c,a} =& \sum_{k=0}^{c} \exp(- A_{k,a}) 
    #         =& B_{c-1,a} + \exp(- A_{c,a}).


    B = util.make_matrix(ntimes, ntimes, 0.0)
    for b in xrange(ntimes):
        B[0][b] = nbranches[0] * time_steps[0] / nrecombs[0] * exp(-A[0][b])
        for c in xrange(1, b):
            B[c][b] = (B[c-1][b] + nbranches[c] * time_steps[c] / nrecombs[c]
                       * exp(-A[c][b]))

    # S_{a,b} &= B_{min(a-1,b-1),a}
    S = util.make_matrix(ntimes, ntimes, 0.0)
    for a in xrange(1, ntimes):
        for b in xrange(1, ntimes):
            S[a][b] = B[min(a-1, b-1)][b]

    # f =\frac{[1 - \exp(- \rho (|T^{n-1}_{i-1}| + s_a))] 
    #       [1 - \exp(- s'_{b-1} k_{b-1} / (2N))]}
    #      {\exp(-\rho |T^{n-1}_{i-1}|) (|T^{n-1}_{i-1}| + s_a) k^C_b}
    # |T^{n-1}_{i-1}| = treelen

    # TODO: fix for case where b=0
    
    time_lookup = util.list2lookup(times)
    transprob = util.make_matrix(len(states), len(states), 0.0)
    for i, (node1, a) in enumerate(states):
        c = time_lookup[tree[node1].age]
        treelen2 = treelen + max(times[a], mintime)
        
        for j, (node2, b) in enumerate(states):
            f = ((1.0 - exp(-rho * treelen2)) /
                 (exp(-rho * treelen) * treelen2 * ncoals[b]))
            if b > 0:
                f *= (1.0 - exp(-time_steps[b-1] * nbranches[b-1]
                                / (2.0 * popsizes[b-1])))
            else:
                # HACK
                f *= 0.0
            if node1 != node2:
                transprob[i][j] = f * S[a][b]
            elif a != b:
                transprob[i][j] = f * (2*S[a][b] - S[c][b])
            else:
                # compute at the end
                pass

        transprob[i][i] = 1.0 - sum(transprob[i])
        for j in xrange(len(states)):
            transprob[i][j] = util.safelog(transprob[i][j])

    return transprob



def calc_transition_probs2_c(tree, states, nlineages, times,
                            time_steps, popsizes, rho, raw=True):
    
    nbranches, nrecombs, ncoals = nlineages

    times_lookup = dict((t, i) for i, t in enumerate(times))
    tree2 = tree.get_tree()
    ptree, nodes, nodelookup = make_ptree(tree2)
    int_states = [[nodelookup[tree2[node]], timei]
                  for node, timei in states]
    nstates = len(int_states)
    ages_index = [times_lookup[tree[node.name].age]
                  for node in nodes]
    #treelen = sum(x.dist for x in tree2)
    treelen = get_treelen(tree, times)
    transmat = new_transition_probs2(
        len(nodes), ages_index, treelen, 
        ((c_int * 2) * nstates)
        (* ((c_int * 2)(n, t) for n, t in int_states)), nstates,
        len(time_steps), times, time_steps,
        nbranches, nrecombs, ncoals, 
        popsizes, rho)

    if raw:
        return transmat
    else:
        transmat2 = [transmat[i][:nstates]
            for i in range(nstates)]
        delete_transition_probs(transmat, nstates)
        return transmat2



def get_deterministic_transitions_debug(states1, states2, times,
                                  tree, last_tree,
                                  recomb_branch, recomb_time,
                                  coal_branch, coal_time):

    # recomb_branch in tree and last_tree
    # coal_branch in last_tree

    def walk_up(node, start, time, ignore=None):
        if (coal_branch == node or coal_branch == start) and coal_time < time:
            # coal occurs under us
            # TODO: make this probabilistic
            ptr = tree2[start].parents[0]
            while len(ptr.children) != 2 or ptr.name == ignore:
                ptr = ptr.parents[0]
            return ptr.name
        else:
            return start

    
    state2_lookup = util.list2lookup(states2)
    last_tree2 = last_tree.copy()
    arglib.remove_single_lineages(last_tree2)
    tree2 = tree.copy()
    arglib.remove_single_lineages(tree2)
    
    
    next_states = []
    for i, state1 in enumerate(states1):
        node1, a = state1
        time = times[a]
        
        if (node1, a) == (coal_branch, coal_time):
            # not a deterministic case
            next_states.append(-1)
        
        elif node1 != recomb_branch:
            # SPR only removes a subset of descendents, if any
            # trace up from remaining leaf to find correct new state

            node = last_tree2.nodes.get(node1, None)
            if node is None:
                print node1
                treelib.draw_tree_names(last_tree.get_tree(),
                                        minlen=8, maxlen=8)
                raise Exception("unknown node name '%s'" % node1)

            
            if node.is_leaf():
                # SPR can't disrupt leaf branch
                node2 = walk_up(node1, node1, a)
                next_states.append(state2_lookup[(node2, a)])

            else:
                child1 = node.children[0]
                while len(child1.children) == 1:
                    child1 = child1.children[0]

                child2 = node.children[1]
                while len(child2.children) == 1:
                    child2 = child2.children[0]
                
                if recomb_branch == child1.name:
                    # right child is not disrupted
                    node2 = walk_up(node1, child2.name, a, node1)
                    next_states.append(state2_lookup[(node2, a)])

                elif recomb_branch == child2.name:
                    # left child is not disrupted
                    node2 = walk_up(node1, child1.name, a, node1)
                    next_states.append(state2_lookup[(node2, a)])

                else:
                    # node is not disrupted
                    node2 = walk_up(node1, node1, a)
                    next_states.append(state2_lookup[(node2, a)])
                  
                
        else:
            # SPR is on same branch as new chromosome
            if recomb_time >= a:
                # we move with SPR subtree
                # TODO: we could probabilistically have subtree move
                # out from underneath.
                next_states.append(state2_lookup[(recomb_branch, a)])

            else:
                # SPR should not be able to coal back onto same branch
                # this would be a self cycle
                if coal_branch == node1:
                    print (recomb_branch, recomb_time), \
                          (coal_branch, coal_time)
                    treelib.draw_tree_names(last_tree.get_tree(),
                                            minlen=8, maxlen=8)
                    treelib.draw_tree_names(tree.get_tree(),
                                            minlen=8, maxlen=8)

                    print "path1"
                    ptr = last_tree[recomb_branch]
                    ptr = ptr.parents[0]
                    while len(ptr.children) == 1:
                        print ptr.name, ptr.event
                        ptr = ptr.parents[0]

                    print "path2"
                    ptr = tree[recomb_branch]
                    ptr = ptr.parents[0]
                    while len(ptr.children) == 1:
                        print ptr.name, ptr.event
                        ptr = ptr.parents[0]
                    
                    assert False

                
                # SPR subtree moves out from underneath us
                # therefore therefore the new chromosome coalesces with
                # the branch above the subtree

                # search up for parent
                recomb = last_tree2[recomb_branch]
                parent = recomb.parents[0]
                b = times.index(parent.age)

                # find other child
                c = parent.children
                other = c[0] if c[1] == recomb else c[1]

                # find new state in tree
                if other.name == coal_branch:
                    next_state = (tree2[other.name].parents[0].name, b)
                else:
                    next_state = (other.name, b)
                
                next_states.append(state2_lookup[next_state])

                # search up for parent
                ptr = last_tree[recomb_branch]
                ptr = ptr.parents[0]
                while len(ptr.children) == 1:
                    ptr = ptr.parents[0]
                b = times.index(ptr.age)

                if ptr.name not in tree:
                    # we are above root
                    assert ptr.age >= tree.root.age
                    ptr = tree.root
                else:
                    ptr = tree[ptr.name]
                    
                # walk down for next coal node
                while len(ptr.children) == 1:
                    ptr = ptr.children[0]

                next_state2 = (ptr.name, b)
                assert next_state2 == next_state

    return next_states



def get_nlineages_recomb_coal2(tree, times):
    """
    Count the number of lineages at each time point that can coal and recomb
    """

    # TODO: add recomb points at end of branches too.
    
    nlineages = [0 for i in times]
    nlineages_recomb = [0 for i in times]
    nlineages_coal = [0 for i in times]

    for name, timei in iter_coal_states(tree, times):
        node = tree[name]

        # find parent node
        if node.parents:
            parent = node.parents[0]
            while len(parent.children) == 1:
                parent = parent.parents[0]
        else:
            parent = None

        # count who passes through this time segment
        if not parent or times[timei] < parent.age:
            nlineages[timei-1] += 1

        # count as recomb unless it is last time point on branch
        if not parent or times[timei] < parent.age:
            nlineages_recomb[timei] += 1

        # count as coal point
        nlineages_coal[timei] += 1
    nlineages[-1] = 1
    
    return nlineages, nlineages_recomb, nlineages_coal



def backward_algorithm2(model, n, verbose=False):

    probs = []

    # calc last position
    nstates = model.get_num_states(n-1)
    for i in xrange(n):
        probs.append(None)
    probs[n-1] = [model.prob_prior(n-1, j) + model.prob_emission(n-1, j)
                  for j in xrange(nstates)]
    
    if n > 20:
        step = (n // 20)
    else:
        step = 1
    
    # loop through positions
    nstates1 = nstates
    i = n-2
    next_print = n-step
    if verbose:
        util.tic("backward")
    while i > -1:
        if verbose and i < next_print:
            next_print -= step
            print " backward iter=%d/%d" % (i+1, n)

        # do first position manually
        nstates1 = model.get_num_states(i)
        nstates2 = model.get_num_states(i+1)
        col2 = probs[i+1]

        model.check_local_tree(i+1)
        if i+1 == model.local_block[0] and model.transmat_switch:
            trans = model.transmat_switch
        else:
            trans = model.transmat


        # find total transition and emission
        col1 = []
        emit = [model.prob_emission(i+1, k) for k in xrange(nstates2)]
        for j in xrange(nstates1):
            tot = -util.INF
            for k in xrange(nstates2):
                p = col2[k] + emit[k] + trans[j][k]
                tot = logadd(tot, p)
            col1.append(tot)
        probs[i] = col1
        i -= 1
        if i <= -1:
            break

        # do rest of block quickly
        space = model.get_state_space(i)
        block = model.get_local_block(space)
        blocklen = i+1 - block[0]
        if i < block[1]-1 and blocklen > 4:
            nstates = model.get_num_states(i)

            # setup tree and states
            tree = model.arg.get_marginal_tree(i-.5)
            tree2 = tree.get_tree()
            ptree, nodes, nodelookup = make_ptree(tree2)
            int_states = [[nodelookup[tree2[node]], timei]
                          for node, timei in model.states[i]]
            ages = [tree[node.name].age for node in nodes]
            seqs = [model.seqs[node.name][block[0]:i+2]
                    for node in nodes if node.is_leaf()]
            seqs.append(model.seqs[model.new_name][block[0]:i+2])
            seqlen = blocklen + 1
            
            emit = new_emissions(
                ((c_int * 2) * nstates)
                (* ((c_int * 2)(n, t) for n, t in int_states)), nstates, 
                ptree, len(ptree), ages,
                (c_char_p * len(seqs))(*seqs), len(seqs), seqlen,
                model.times, len(model.times), model.mu)
            
            trans = c_matrix(c_double,
                             [[model.prob_transition(i, j, i+1, k)
                               for k in xrange(nstates)]
                              for j in xrange(nstates)])
            bw = [[0.0 for k in xrange(nstates)]
                  for pos in xrange(block[0], i+1)]
            bw.append(probs[i+1])
            
            backward_alg(blocklen+1, nstates, trans, emit, bw)
            for j in xrange(blocklen):
                probs[block[0]+j] = bw[j][:nstates]
            i = block[0] - 1

    if verbose:
        util.toc()

    return probs



'''
def get_nlineages(tree, times):
    """Count the number of lineages in each time segment"""
    nlineages = [0 for i in times]
    for name, i in iter_coal_states(tree, times):
        node = tree[name]
        if node.parents:
            parent = node.parents[0]
            while len(parent.children) == 1:
                parent = parent.parents[0]
        if not node.parents or times[i] < parent.age:
            nlineages[i-1] += 1
    nlineages[-1] = 1
    return nlineages
'''



'''
def add_arg_thread3(arg, new_name, thread, recombs):

    # XXX: cycles are being made

    given = set(arg)
    
    def walk_up(node, pos, time):
        parent = arg.get_local_parent(node, pos-.5)
        while parent and parent.age < time:
            node = parent
            parent = arg.get_local_parent(node, pos-.5)
        return node, parent

    def walk_up2(node, pos, time):
        parent = arg.get_local_parent(node, pos-.5)
        while parent and parent.age <= time:
            node = parent
            #print "  ...", node.name, node.event, node.pos, node.age, node in given
            parent = arg.get_local_parent(node, pos-.5)
        return node, parent

    def add_parent(node, parent, pos):
        if node.event == "recomb":
            if pos < node.pos:
                assert node.parents[0] is None
                node.parents[0] = parent
            else:
                assert node.parents[1] is None
                node.parents[1] = parent
        else:
            node.parents.append(parent)


    # get postorder index of nodes in arg
    node_order = dict((node.name, i) for i, node in enumerate(arg.postorder()))

    # collect all events in order of age
    events = []

    # get recombination events
    # we subtract 1 from pos because of the difference in notation
    for pos, recomb_node, recomb_time in recombs:
        events.append(("recomb", recomb_time, recomb_node, pos-1))

    # group coal points
    """
    coals = [[]]
    last = None
    j = 0
    for pos, coal_point in enumerate(thread):
        if j < len(recombs) and pos == recombs[j][0]:
            assert coal_point[1] >= recombs[j][2]
            coals.append([])
            last = None
            j += 1
        if coal_point != last:
            coals[-1].append(("coal", coal_point[1], coal_point[0], pos))
            last = coal_point
    """

    last = None
    coal_orders = {}
    coal_order = 0
    j = 0
    arg_recombs  = set(x.pos for x in arg if x.event == "recomb")
    for pos, coal_point in enumerate(thread):
        if j < len(recombs) and pos == recombs[j][0]-1:
            last = None
        if j < len(recombs) and pos == recombs[j][0]:
            assert coal_point[1] >= recombs[j][2]
            last = None
            j += 1
        if pos in arg_recombs or (pos-1) in arg_recombs:
            last = None
        if coal_point != last:
            events.append(
                ("coal", coal_point[1], coal_point[0], pos, coal_order))
            last = coal_point

    print "UNSORTED"
    for event in events:
        print event
    print


    events.sort(key=lambda x: (x[1], x[0] == "coal", x[3]))

    print "ncoals", util.count(lambda x: x[0] == "coal", events)
    print "nrecomb", util.count(lambda x: x[0] == "recomb", events)

    # start new chromsome lineage
    leaf = arg.new_node(name=new_name, age=0.0, event="gene")
    leaf.data["ancestral"] = [(arg.start, arg.end)]
    nstubs = 1
    
    
    # process events
    for event in events:
        print "event", event
        print "nstubs", nstubs
        
        if event[0] == "recomb":
            # recomb event
            
            tag, recomb_time, recomb_name, pos = event

            # find node that is recombining
            node = arg[recomb_name]
            node, parent = walk_up2(node, pos, recomb_time)

            # create recomb node
            recomb = arg.new_node(age=recomb_time, event="recomb",
                                  pos=pos, children=[node],
                                  parents=[None, None])
            
            if parent:
                print "* break"
                # break branch
                node.parents[node.parents.index(parent)] = recomb
                parent.children[parent.children.index(node)] = recomb

                if recomb_name == new_name:
                    recomb.parents[1] = parent
                else:
                    recomb.parents[0] = parent
            else:
                # off local tree
                add_parent(node, recomb, pos)

            nstubs += 1
                

        elif event[0] == "coal":
            # coal event
            # NOTE: coal pos needs -.5 when using add_parent()

            tag, coal_time, coal_name, pos, o = event
            
            # find new_node lineage
            node1, parent1 = walk_up2(leaf, pos, coal_time)
            
            # find node that new_node coals with
            node2 = arg[coal_name]
            node2, parent2 = walk_up2(node2, pos, coal_time)

            if node1 == node2:
                print "SKIP"
                continue

            if parent1 and parent2:
                print "SKIP2"
                continue


            # create coal node
            coal = arg.new_node(age=coal_time, event="coal",
                                children=[node1, node2])

            
            if parent1:
                print "* break1"
                # break branch
                node1.parents[node1.parents.index(parent1)] = coal
                parent1.children[parent1.children.index(node1)] = coal
                coal.parents.append(parent1)
            else:
                # coal is off ARG
                add_parent(node1, coal, pos-.5)

            # is there a chance of the parent changing because of break with
            # node1 and parent1
            if parent2:
                print "* break2"
                # break branch
                node2.parents[node2.parents.index(parent2)] = coal
                parent2.children[parent2.children.index(node2)] = coal
                coal.parents.append(parent2)
            else:
                # coal is off ARG
                add_parent(node2, coal, pos-.5)

            nstubs -=1 

    for node in arg:
        for i, parent in enumerate(node.parents):
            if parent is None:
                print "NULL", node, i, node.event, node.pos, node.age
'''            

    


'''
def add_arg_thread2(arg, new_name, thread, recombs):

    # ensure ancestral sequence is computed

    def find_thread_stub(thread_stubs, pos):
        lineage = None
        for stub in thread_stubs:
            if stub[0] < pos-.5 < stub[1]:
                lineage = stub
                break
        assert lineage is not None
        return lineage

    def walk_up(node, pos, time):
        parent = arg.get_local_parent(node, pos-.5)
        while parent and parent.age <= time:
            node = parent
            parent = arg.get_local_parent(node, pos-.5)
        return node, parent

    def add_parent(node, parent, pos):
        if node.event == "recomb":
            if pos < node.pos:
                assert node.parents[0] == None
                node.parents[0] = parent
            else:
                assert node.parents[1] == None
                node.parents[1] = parent
        else:
            node.parents.append(parent)


    # get postorder index of nodes in arg
    node_order = dict((node.name, i) for i, node in enumerate(arg.postorder()))

    # collect all events in order of age
    events = []

    # get recombination events
    for pos, recomb_node, recomb_time in recombs:
        events.append(("recomb", recomb_time, recomb_node, pos))

    # group coal points
    coals = [[]]
    last = None
    j = 0
    for pos, coal_point in enumerate(thread):
        if j < len(recombs) and pos == recombs[j][0] + 1:
            coals.append([])
            last = None
            j += 1
        if coal_point != last:
            if pos > 0:
                assert coal_point[1] >= recombs[j-1][2], (pos, j-1)
            coals[-1].append(("coal", coal_point[1], coal_point[0], pos))
            last = coal_point

    # determine coal event per group
    for coal_group in coals:
        # find min age coals
        minage = min(x[1] for x in coal_group)
        events.append(max((x for x in coal_group if x[1] == minage),
                          key=lambda x: node_order[x[2]]))

    print coals

    events.sort(key=lambda x: (x[1], x[3]))

    print "ncoals", util.count(lambda x: x[0] == "coal", events)
    print "nrecomb", util.count(lambda x: x[0] == "recomb", events)

    # start new chromsome lineage
    leaf = arg.new_node(name=new_name, age=0.0, event="gene")
    leaf.data["ancestral"] = [(arg.start, arg.end)]

    thread_stubs = set([(-1, arg.end, 0, leaf)])
    thread_coal = set()
    other_stubs = 0
    
    # process events
    for event in events:
        print
        print "event", event
        print "stubs", thread_stubs
        print "coaled", thread_coal
        
        if event[0] == "recomb":
            # recomb event
            
            tag, recomb_time, recomb_name, pos = event
            if recomb_name == new_name:
                # recomb in new thread lineage

                # find lineage to recombine
                lineage = find_thread_stub(thread_stubs, pos)
                thread_stubs.remove(lineage)
                start, end, side, node = lineage

                # create recomb node
                recomb = arg.new_node(age=recomb_time, event="recomb",
                                      pos=pos, children=[node],
                                      parents=[None, None])
                add_parent(node, recomb, pos)
                    
                # create two new lineages
                thread_stubs.add((start, pos, 0, recomb))
                thread_stubs.add((pos, end, 1, recomb))
            else:
                # recomb in G_{n-1}

                # find node and region in G_{n-1} that is recombining
                node = arg[recomb_name]
                node, parent = walk_up(node, pos, recomb_time)
                
                # create recomb node
                recomb = arg.new_node(age=recomb_time, event="recomb",
                                      pos=pos, children=[node],
                                      parents=[None, None])
                
                if parent:
                    # break branch
                    node.parents[node.parents.index(parent)] = recomb
                    parent.children[parent.children.index(node)] = recomb
                    recomb.parents[0] = parent
                    other_stubs += 1
                else:
                    # off local tree
                    add_parent(node, recomb, pos)
                    other_stubs += 1
                        

        elif event[0] == "coal":
            # coal event
            # NOTE: coal pos needs -.5 when using add_parent()

            tag, coal_time, coal_name, pos = event

            # find stub in thread that is coalescing
            lineage = find_thread_stub(thread_stubs, pos)
            start1, end1, side1, node1 = lineage

            # find node in G_{n-1} to coalesce with
            node2 = arg[coal_name]
            node2, parent2 = walk_up(node2, pos, coal_time)
            
            if (start1, end1) in thread_coal:
                print "thread is coaled"
                # we have already coalesced node1
                # see if additional branches are needed
                if parent2 is not None:       
                    print node2.name, parent2.name
                    treelib.draw_tree_names(
                        arg.get_marginal_tree(pos).get_tree(),
                        minlen=5, maxlen=5)
                    assert False
                

                # walk up to find coal point
                node1, parent1 = walk_up(node1, pos, coal_time)
                
                # create coal node
                coal = arg.new_node(age=coal_time, event="coal",
                                    children=[node1, node2])
                add_parent(node2, coal, pos-.5)
                other_stubs -=1
                
                if parent1:
                    # break branch
                    node1.parents[node1.parents.index(parent1)] = coal
                    parent1.children[parent1.children.index(node1)] = coal
                    coal.parents.append(parent1)
                else:
                    # coal is off ARG
                    add_parent(node1, coal, pos-.5)
                
            else:
                # we have not coalesced node1 yet
                print "thread is not coaled"
                
                # create coal node
                coal = arg.new_node(age=coal_time, event="coal",
                                    children=[node1, node2])
                thread_stubs.remove(lineage)
                thread_stubs.add((start1, end1, 0, coal))

                if node1.event == "recomb":
                    node1.parents[side1] = coal
                else:
                    node1.parents.append(coal)

                if parent2:
                    # break branch
                    node2.parents[node2.parents.index(parent2)] = coal
                    parent2.children[parent2.children.index(node2)] = coal
                    coal.parents.append(parent2)
                    thread_coal.add((start1, end1))
                else:
                    # coal is off ARG
                    add_parent(node2, coal, pos-.5)
                    other_stubs -= 1

            
        print "nstubs", len(thread_stubs) - len(thread_coal), other_stubs
'''                
            
'''

def add_arg_thread3(arg, new_name, thread, recombs, arg3=None):


    def is_local_coal(arg, node, pos, local):
        return (len(node.children) == 2 and
                node.children[0] in local and
                arg.get_local_parent(node.children[0], pos-.5) == node and
                node.children[1] in local and
                arg.get_local_parent(node.children[1], pos-.5) == node and
                node.children[0] != node.children[1])



    def walk_up(arg, leaves, time, pos, ignore=None):

        order = dict((node, i) for i, node in enumerate(
            arg.postorder_marginal_tree(pos-.5)))
        local = set(order.keys())
        if ignore is not None and ignore in arg:
            ptr = arg[ignore]
            local.remove(ptr)
            ptr = arg.get_local_parent(ptr, pos-.5)
            
            while ptr:
                if (len(ptr.children) == 2 and
                    ((ptr.children[0] in local and
                      arg.get_local_parent(ptr.children[0], pos-.5) == ptr) or
                     (ptr.children[1] in local and
                      arg.get_local_parent(ptr.children[1], pos-.5) == ptr))):
                    print "halt", ptr
                    break
                print "remove", ptr
                local.remove(ptr)
                ptr = arg.get_local_parent(ptr, pos-.5)
        print local

        queue = [(order[arg[x]], arg[x]) for x in leaves]
        seen = set(x[1] for x in queue)
        heapq.heapify(queue)

        while len(queue) > 1:
            i, node = heapq.heappop(queue)
            parent = arg.get_local_parent(node, pos-.5)
            if parent and parent not in seen:
                seen.add(parent)
                heapq.heappush(queue, (order[parent], parent))
        node = queue[0][1]
        parent = arg.get_local_parent(node, pos-.5)

        
        while parent and parent.age <= time:
            if is_local_coal(arg, parent, pos, local):
                print "stop", node.age, parent.age, parent.children
                break            
            node = parent
            parent = arg.get_local_parent(node, pos-.5)

        if parent:
            if parent.age < time:
                print (leaves, parent.age, time)
                tree = arg.get_marginal_tree(pos-.5).get_tree()
                tree.write()
                treelib.draw_tree_names(tree, maxlen=8, minlen=8)
                assert False

        return node


    def add_node(arg, node, time, pos, event):

        node2 = arg.new_node(event=event, age=time, children=[node], pos=pos)
        if event == "coal":
            node2.pos = 0

        parent = arg.get_local_parent(node, pos-.5)
        if parent:
            node.parents[node.parents.index(parent)] = node2
            parent.children[parent.children.index(node)] = node2
            node2.parents.append(parent)
        else:
            node.parents.append(node2)

        return node2


    arg_recomb = dict((x.pos, x) for x in iter_visible_recombs(arg))
    recomb_clades = [
        (pos-1, None) + get_clade_point(arg, rnode, rtime, pos-1)
        for pos, rnode, rtime in recombs] + [
        (node.pos, node.name) +
        get_clade_point(arg, node.name, node.age, node.pos)
        for node in iter_visible_recombs(arg)]
    recomb_clades.sort()

    # make initial tree
    arg2 = arg.get_marginal_tree(-1)
    arglib.remove_single_lineages(arg2)

    start = get_clade_point(arg, thread[0][0], thread[0][1], 0)
    node = walk_up(arg2, start[0], start[1], -1)
    node2 = add_node(arg2, node, start[1], -1, "coal")
    leaf = arg2.new_node(name=new_name, event="gene", age=0)
    leaf.parents.append(node2)
    node2.children.append(leaf)

    arg2.get_marginal_tree(-.5).get_tree().write()


    # add each recomb and re-coal
    for rpos, rname, rleaves, rtime in recomb_clades:
        print "------------------------------------------"

        for node in arg2:
            if node.event == "recomb":
                assert len(node.parents) == 2, node
        
        if rpos in arg_recomb:
            # find re-coal for existing recomb

            node = arg_recomb[rpos]
            local1 = set(arg.postorder_marginal_tree(rpos-.5))
            local2 = set(arg.postorder_marginal_tree(rpos+.5))
            last = node
            node = arg.get_local_parent(node, rpos+.5)
            while (node not in local1 and
                   not is_local_coal(arg, node, rpos+1, local2)):
                last = node
                node = arg.get_local_parent(node, rpos+.5)
            c = node.children
            child = c[0] if c[1] == last else c[1]

            print ">>", node, c
            treelib.draw_tree_names(
                arg.get_marginal_tree(rpos-.5).get_tree(),
                maxlen=8, minlen=8)
            treelib.draw_tree_names(
                arg.get_marginal_tree(rpos+.5).get_tree(),
                maxlen=8, minlen=8)
                
            cleaves, ctime = get_clade_point(
                arg, child.name, node.age, rpos-.5)


            # determine if this is mediated re-coal
            #  1. previous coal-state == re-coal point
            #  2. this coal-state == is above recomb node
            leaves1, time1 = get_clade_point(arg, thread[rpos][0],
                                             thread[rpos][1], rpos-.5)
            leaves2, time2 = get_clade_point(arg, thread[rpos+1][0],
                                             thread[rpos+1][1], rpos+1-.5)

            if (set(leaves1) == set(cleaves) and time1 == ctime and
                set(leaves2) == set(rleaves)):
                # mediated
                print "mediated", (leaves1, time1), (leaves2, time2)
                print "  ", (rleaves, rtime), (cleaves, ctime)
                
                cleaves = [new_name]
                ctime = thread[rpos+1][1]

                node1 = walk_up(arg2, rleaves, rtime, rpos+1, new_name)
                node2 = walk_up(arg2, cleaves, ctime, rpos+1)
            else:
                node1 = walk_up(arg2, rleaves, rtime, rpos+1, new_name)
                node2 = walk_up(arg2, cleaves, ctime, rpos+1, new_name)
    
        else:
            # find re-coal for new recomb
            
            assert rtime <= thread[rpos][1], (rtime, thread[rpos][1])
            
            if rleaves == [new_name]:
                # recomb on new branch, coal given thread
                cleaves, ctime = get_clade_point(
                    arg, thread[rpos+1][0], thread[rpos+1][1], rpos+.5)
                assert ctime >= rtime, (rtime, ctime)
                
                node1 = walk_up(arg2, rleaves, rtime, rpos+1)
                node2 = walk_up(arg2, cleaves, ctime, rpos+1, new_name)
                
            else:
                # recomb in ARG, coal on new branch
                cleaves = [new_name]
                ctime = thread[rpos+1][1]
                assert ctime >= rtime, (rtime, ctime)
                
                node1 = walk_up(arg2, rleaves, rtime, rpos+1, new_name)
                node2 = walk_up(arg2, cleaves, ctime, rpos+1, node1.name)



        print "add", rpos, rpos in arg_recomb, (rleaves, rtime), (cleaves, ctime)

        if arg3:
            print "arg3"
            tree = arg3.get_marginal_tree(rpos+.5).get_tree()
            treelib.remove_single_children(tree)
            phylo.hash_order_tree(tree)
            tree.write()

        print "arg"
        tree = arg.get_marginal_tree(rpos+.5).get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        tree.write()

        print "arg (i-1)"
        tree = arg.get_marginal_tree(rpos-.5).get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        tree.write()

        #if set(rleaves) == set(cleaves):
        #    print "skip", rleaves
        #    continue

        assert node1.parents
        assert rtime <= ctime

        recomb = add_node(arg2, node1, rtime, rpos, "recomb")
        if node1 == node2:
            node2 = recomb
        coal = add_node(arg2, node2, ctime, rpos, "coal")

        recomb.parents.append(coal)
        coal.children.append(recomb)

        print "arg2"
        tree = arg2.get_marginal_tree(rpos+.5).get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        tree.write()

        print "arg2 (n-1)"
        tree = arg2.get_marginal_tree(rpos-.5).get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        tree.write()


        print "  r", recomb, recomb.children, recomb.parents
        print "  c", coal, coal.children, coal.parents



        arglib.assert_arg(arg2)

        node, time = get_coal_point(arg2, arg2[new_name], rpos+1)
        assert time == thread[rpos+1][1], (time, thread[rpos+1][1])

    
    
    return arg2
    
'''


'''
def get_deterministic_transition(state1, states2, times, tree, last_tree,
                                 recomb_branch, recomb_time,
                                 coal_branch, coal_time):

    # recomb_branch in tree
    # coal_branch in last_tree
    
    node1, a = state1

    # get leaves under node1
    leaves1 = set(last_tree.leaves(last_tree[node1]))

    # get leaves under recomb_node
    recomb_leaves = set(last_tree.leaves(last_tree[recomb_branch]))

    remain = leaves1 - recomb_leaves

    def find_state(node, time):
        b = util.INF
        state2 = None
        
        while len(node.children) == 1:
            node = node.children[0]
        
        for j, (n, t) in enumerate(states2):
            if node.name == n and time <= t < b:
                b = t
                state2 = j
        assert state2 is not None, ((node, time), states2)
        return state2
                

    def trace_up(node, time):
        last = node
        while node.age <= times[time]:
            if len(node.children) != 1:
                last = node
            if not node.parents:
                break
            node = node.parents[0]
        return last
    
    
    if len(remain) > 0:
        # SPR only removes a subset of descendents
        # trace up from remaining leaf to find correct new state
        ptr = tree[remain.pop().name]
        node = trace_up(ptr, a)
        return find_state(node, a)

    else:
        # SPR is on same branch as new chromosome
        if recomb_time >= a:
            # we move with SPR subtree
            ptr = tree[recomb_leaves.pop().name]
            node = trace_up(ptr, a)
            return find_state(node, a)

        elif coal_time <= a and coal_branch == node1:
            # SPR subtree coals back underneath us
            return find_state(tree[node1], a)
            
        else:
            # SPR subtree moves out from underneath us
            # therefore therefore the new chromosome coalesces with
            # the branch above the subtree

            # search up for parent
            ptr = last_tree[recomb_branch]
            ptr = ptr.parents[0]
            while len(ptr.children) == 1:
                ptr = ptr.parents[0]
            b = times.index(ptr.age)

            # go over to new tree
            if ptr.name not in tree:
                # we are above root
                assert ptr.age >= tree.root.age
                return find_state(tree.root, b)
            ptr = tree[ptr.name]
            
            return find_state(ptr, b)
'''        


'''
def get_deterministic_transition2(state1, states2, times, tree, last_tree,
                                 recomb_branch, recomb_time,
                                 coal_branch, coal_time):

    # recomb_branch in tree
    # coal_branch in last_tree

    node1, a = state1

    if recomb_branch == node1:
        if recomb_time >= a:
            # state unaffected
            return states2.index(state1)
        else:
            # search up for parent
            ptr = last_tree[recomb_branch]
            last_ptr = ptr
            ptr = ptr.parents[0]
            while len(ptr.children) == 1:
                last_ptr = ptr
                ptr = ptr.parents[0]
            b = times.index(ptr.age)

            # go over to new tree
            if ptr.name not in tree:
                # we are above root
                assert ptr.age >= tree.root.age
                print tree.root.name, b
                # may have to go up coal point is not available
                #while (tree.root.name, b) not in states2:
                #    b += 1
                return states2.index((tree.root.name, b))
            ptr = tree[ptr.name]

            # go down to beginning of branch
            while len(ptr.children) == 1:
                ptr = ptr.children[0]

            # may have to go up coal point is not available
            while (ptr.name, b) not in states2:
                b += 1
            print ptr.name, b
            return states2.index((ptr.name, b))

    elif coal_branch == node1 and coal_time < a:
        # move to coal branch
        ptr = tree[coal_branch]
        last = ptr
        print last.name
        while ptr.age <= times[a] and ptr.parents:
            if len(ptr.children) != 1:
                last = ptr
            ptr = ptr.parents[0]
        print last.name, a
        # may have to go up coal point is not available
        #while (last.name, a) not in states2:
        #   a += 1
        return states2.index((last.name, a))

    elif coal_branch == node1 and coal_time == a:
        raise Exception("transition is not deterministic")

    else:
        # either unaffected or moved to next node with 0 or 2 children
        if node1 in tree:
            ptr = tree[node1]
        else:
            ptr = tree.root
        while len(ptr.children) == 1:
            ptr = ptr.children[0]
        return states2.index((ptr.name, a))

'''



'''
def sample_recombinations(model, path, use_times=True):
    
    r = 0
    # assumes that recomb_pos starts with -1 and ends with arg.end
    arg_recomb = model.recomb_pos
    tree = model.arg.get_marginal_tree(0)
    treelen = sum(x.get_dist() for x in tree)
    new_node = model.new_name
    
    for pos, state in enumerate(path):
        node, timei = model.states[pos][state]
        #print r, arg_recomb[r]
        
        # update local tree if needed
        while r < len(arg_recomb) and arg_recomb[r] < pos:
            r += 1
            model.check_local_tree(pos)
            tree = model.local_tree
            treelen = sum(x.get_dist() for x in tree)
            nbranches = model.nlineages[0]
            nrecombs = model.nlineages[1]
            ncoals =  model.nlineages[2]
            A = calc_A_matrix(model.time_steps, nbranches, model.popsizes)

        if pos == 0 or arg_recomb[r-1] == pos - 1:
            # previous arg recomb is right behind us, sample no recomb
            continue

        # get information about pos-1
        # since their no recomb in G_{n-1}, last_tree == tree
        last_state = path[pos-1]
        last_node, last_timei = model.states[pos-1][last_state]
        last_tree = tree
        last_treelen = treelen

        blen = model.times[last_timei]
        last_treelen2 = last_treelen + blen
        if node == last_tree.root.name:
            last_treelen2 += blen - last_tree.root.age

        if state == last_state:
            # state is the same, there is a chance of no recomb
            p = exp(-model.rho * (last_treelen2 - last_treelen))

            if random.random() < p:
                # sample no recombination
                continue

        # there must be a recombination
        # either because state changed or we choose to recombine

        if node == last_node:
            if timei == last_timei:
                # y = v, k in [0, min(timei, last_timei))
                # y = node, k in Sr(node)
                # if node.parent.age == model.times[timei],
                #   y = sis(last_tree, node.name), k in Sr(y)
                node_timei = model.times.index(tree[node].age)
                recombs = [(new_node, k) for k in
                           range(0, min(timei, last_timei))] + \
                          [(node, k) for k in
                           range(node_timei, min(timei, last_timei))]

                # TODO: add back
                #if last_tree[node].parents[0].age == model.times[timei]:
                #    y = None # TODO: sis(node)
                #    recombs += [(y, k) for k in
                #                range(model.times.index(y.age),
                #                      min(timei, last_timei))]
                
            else:
                # y = v, k in [0, min(timei, last_timei))
                # y = node, k in Sr(node)
                node_timei = model.times.index(tree[node].age)
                recombs = [(new_node, k) for k in
                           range(0, min(timei, last_timei))] + \
                          [(node, k) for k in
                           range(node_timei, min(timei, last_timei))]
            
        else:
            # y = v, k in [0, min(timei, last_timei))
            recombs = [(new_node, k) for k in range(0, min(timei, last_timei))]

        if len(recombs) == 0:
            print "SKIP"
            continue
        
        j = timei
        probs = []
        for recomb in recombs:
            probs.append((nbranches[k] + 1) * model.time_steps[k] /
                         (ncoals[j] * (nrecombs[k] + 1) * last_treelen2) *
                         (1.0 - exp(-model.time_steps[j-1] * nbranches[j-1] /
                                    (2.0 * model.popsizes[j-1]))) *
                         (1.0 - exp(-model.rho * last_treelen2)) *
                         exp(-A[k][j]))
        recomb_node, recomb_time = recombs[stats.sample(probs)]
        if use_times:
            recomb_time = model.times[recomb_time]
        yield (pos, recomb_node, recomb_time)
'''


'''

def topsort(graph):

    openset = set()
    visited = {}
    seen = set()

    # assert graph is valid
    for name, node in graph.iteritems():
        assert node.name == name, name
        for c in node.children:
            assert name in graph[c].parents, name
        for p in node.parents:
            assert name in graph[p].children, name 
        

    for name, node in graph.iteritems():
        visited[name] = 0
        if not node.children:
            openset.add(node)

    while openset:
        n = openset.pop()
        assert n.name not in seen, n.name
        seen.add(n.name)
        yield n.name

        for m_name in n.parents:
            m = graph[m_name]
            visited[m_name] += 1
            #m.children.remove(n.name)
            if visited[m_name] == len(m.children):
                assert m_name not in seen, m_name
                openset.add(m)

    for n in graph.values():
        assert len(n.children) == visited[n.name], \
               (n.name, n.children, visited[n.name], n.name in seen)


def find_cycle(graph, events, times, arg):

    def print_cycle(name, seen):
        print name, seen

    def walk(node, seen):
        print "w", node.name, len(seen)
        if node.name[0] == "event":
            e = events[node.name[1]]
            print "  ", e, times.index(e[1])
        elif node.name[0] == "arg":
            n = arg[node.name[1]]
            t = -1 
            if n.age in times:
                t = times.index(n.age)
            print "  ", n.name, n.event, n.pos, t
            
        for p in node.parents:
            if p in seen:
                print_cycle(p, seen)
                assert False
            seen.append(p)
            walk(graph[p], seen)
            seen.pop()

    for name, node in graph.items():
        if len(node.children) == 0:
            walk(node, [node.name])
    


def is_same_coal_event(state1, state2, arg, pos):

    def walk_up(arg, node, pos, node2):
        parent = arg.get_local_parent(node, pos-.5)
        while parent:
            if parent == node2:
                return True
            parent = arg.get_local_parent(parent, pos-.5)
        return False

    node1, a = state1
    node2, b = state2
    
    if a != b:
        return False, -1

    if node1 == node2:
        return True, 0

    print "more", pos

    if walk_up(arg, arg[node1], pos - 1, arg[node2]):
        return True, -1
    if walk_up(arg, arg[node2], pos, arg[node1]):
        return True, -1

    print "mediated", pos

    return False, -1



    tree = arg.get_marginal_tree(pos - .5)
    last_tree = arg.get_marginal_tree(pos - 1.5)


    (recomb_branch, recomb_time), (coal_branch, coal_time) = \
        find_recomb_coal(tree, last_tree, pos=pos)
    # recomb_branch in tree
    # coal_branch in last_tree
    
    # get leaves under recomb_node
    recomb_leaves = set(last_tree.leaves(last_tree[recomb_branch]))    
    leaves1 = set(last_tree.leaves(last_tree[node1]))
    remain = leaves1 - recomb_leaves

    if len(remain) > 0:
        # SPR only removes a subset of descendents
        return True

    else:
        # SPR is on same branch as new chromosome
        if recomb_time > a:
            # we move with SPR subtree
            return True

        elif coal_time < a and coal_branch == node1:
            # SPR subtree coals back underneath us
            return True

        else:
            # SPR subtree moves out from underneath us
            # therefore therefore the new chromosome coalesces with
            # the branch above the subtree
            return False



def add_arg_thread(arg, new_name, thread, recombs):
    
    def walk_up2(node, pos, time, local):
        assert node.name in local
        parent = arg.get_local_parent(node, pos-.5)
        while parent:
            if parent.name in local:
                break
            parent = arg.get_local_parent(parent, pos-.5)
        return node, parent

    def walk_up3(node, pos, time, event_name):
        parent = arg.get_local_parent(node, pos-.5)
        stop = set(graph[event_name].parents)
        #print stop

        while parent:
            parent_name = new_nodes.get(parent, ("arg", parent.name))
            #print parent_name, parent.age
            if parent.age > time or parent_name in stop:
                assert time >= node.age
                break
            node = parent
            parent = arg.get_local_parent(node, pos-.5)
        return node, parent


    def add_parent(node, parent, pos):
        if node.event == "recomb":
            if pos < node.pos:
                assert node.parents[0] is None
                node.parents[0] = parent
            else:
                assert node.parents[1] is None
                node.parents[1] = parent
        else:
            node.parents.append(parent)

    class Node (object):
        def __init__(self, name):
            self.name = name
            self.parents = []
            self.children = []


    # collect all events in order of age
    events = []

    # get recombination events
    # we subtract 1 from pos because of the difference in notation
    for pos, recomb_node, recomb_time in recombs:
        events.append(("recomb", recomb_time, recomb_node, pos-1))

    last = None
    j = 0
    arg_recombs = set(x.pos for x in iter_visible_recombs(arg))
    for pos, coal_point in enumerate(thread):
        if j < len(recombs) and pos == recombs[j][0]-1:
            last = None
        if j < len(recombs) and pos == recombs[j][0]:
            assert coal_point[1] >= recombs[j][2]
            last = None
            j += 1
        if pos in arg_recombs or (pos-1) in arg_recombs:
            last = None
        if coal_point != last:
            events.append(
                ("coal", coal_point[1], coal_point[0], pos))
            last = coal_point


    # sort by position
    events.sort(key=lambda x: (x[3], x[0] == "recomb"))


    # topologically sort events using a graph
    graph = {}

    # add time points
    times = sorted(set(e[1] for e in events))
    for i, time in enumerate(times):
        name = ("time", i)
        n = graph[name] = Node(name)
        if i > 0:
            n.children.append(("time", i-1))
            graph[("time", i-1)].parents.append(name)

    # add arg nodes to sorting graph
    for node in arg:
        name = ("arg", node.name)
        n = graph[name] = Node(name)
        for p in node.parents:
            n.parents.append(("arg", p.name))
        for c in node.children:
            n.children.append(("arg", c.name))

    # add events
    coals = {}
    for eventi, event in enumerate(events):
        tag, event_time, node_name, pos = event
        event_name = ("event", eventi)
        n = graph[event_name] = Node(event_name)

        timei = times.index(event_time)
        time_name = ("time", timei)
        n.parents.append(time_name)
        graph[time_name].children.append(event_name)

        if timei > 0:
            time_name = ("time", timei-1)
            n.children.append(time_name)
            graph[time_name].parents.append(event_name)

        if tag == "coal":
            coals[pos] = event_name

    equals = []
    for eventi, event in enumerate(events):
        tag, event_time, node_name, pos = event
        event_name = ("event", eventi)
        n = graph[event_name]

        if pos == 0 or pos-1 in arg_recombs:
            tree = arg.get_marginal_tree(pos-.5)
            local = set(x.name for x in tree if len(x.children) != 1)

            # determine if coal is equivalent to previous one
            if pos == 0:
                same, where_r = False, 0
            else:
                state1 = thread[pos-1]
                state2 = thread[pos]
                same, where_r = is_same_coal_event(state1, state2, arg, pos)
            print "same", same, where_r

            if same:
                equals.append((coals[pos-1], coals[pos]))
            if where_r < 0:
                recomb = [x for x in arg if x.pos == pos-1][0]
                name = ("arg", recomb.name)
                n.children.append(name)
                graph[name].parents.append(event_name)
                

        if tag == "recomb": 
            # a recomb has to be lower than its neighboring coals
            n.parents.append(("event", eventi-1))
            n.parents.append(("event", eventi+1))
            graph[("event", eventi-1)].children.append(event_name)
            graph[("event", eventi+1)].children.append(event_name)
            
        if node_name != new_name:
            node, parent = walk_up2(arg[node_name], pos, event_time, local)

            assert node.age <= event_time

            node_name = ("arg", node.name)
            n.children.append(node_name)
            graph[node_name].parents.append(event_name)
            
            if parent:
                assert event_time <= parent.age, (
                    node.age, event_time, parent.age, event)
                parent_name = ("arg", parent.name)
                n.parents.append(parent_name)
                graph[parent_name].children.append(event_name)



    # apply equals to graph
    for event_name1, event_name2 in equals:
        print "e", event_name1, event_name2
        print "  ", events[event_name1[1]]
        print "  ", events[event_name2[1]]
        parents1 = graph[event_name1].parents
        parents2 = graph[event_name2].parents
        children1 = graph[event_name1].children
        children2 = graph[event_name2].children

        for x in parents1:
            if x not in parents2:
                parents2.append(x)
                graph[x].children.append(event_name2)
        for x in parents2:
            if x not in parents1:
                parents1.append(x)
                graph[x].children.append(event_name1)
        for x in children1:
            if x not in children2:
                children2.append(x)
                graph[x].parents.append(event_name2)
        for x in children2:
            if x not in children1:
                children1.append(x)
                graph[x].parents.append(event_name1)
    

    #find_cycle(graph, events, times, arg)

    # get topological sort
    order = list(topsort(graph))
    print
    print order

    # check order
    for i, name in enumerate(order):
        for p in graph[name].parents:
            assert i < order.index(p)
        for c in graph[name].children:
            assert i > order.index(c)
    
    
    # get sorted events
    event_names = [name for name in order if name[0] == "event"]
    events = [events[name[1]] for name in order
              if name[0] == "event"]


    # use stable sort to establish correct time
    print events
    assert events == sorted(events, key=lambda x: x[1])



    print "ncoals", util.count(lambda x: x[0] == "coal", events)
    print "nrecomb", util.count(lambda x: x[0] == "recomb", events)

    # start new chromsome lineage
    leaf = arg.new_node(name=new_name, age=0.0, event="gene")
    leaf.data["ancestral"] = [(arg.start, arg.end)]
    nstubs = 1

    new_nodes = {}
    
    # process eventsr
    for event, event_name in izip(events, event_names):
        print "event", event
        print "nstubs", nstubs
        
        if event[0] == "recomb":
            # recomb event
            
            tag, recomb_time, recomb_name, pos = event

            # find node that is recombining
            node = arg[recomb_name]
            node, parent = walk_up3(node, pos, recomb_time, event_name)

            # create recomb node
            recomb = arg.new_node(age=recomb_time, event="recomb",
                                  pos=pos, children=[node],
                                  parents=[None, None])
            new_nodes[recomb] = event_name
            
            assert node.age <= recomb_time
            
            if parent:
                # break branch
                node.parents[node.parents.index(parent)] = recomb
                parent.children[parent.children.index(node)] = recomb

                if recomb_name == new_name:
                    recomb.parents[1] = parent
                else:
                    recomb.parents[0] = parent
            else:
                # off local tree
                add_parent(node, recomb, pos)

            nstubs += 1
                

        elif event[0] == "coal":
            # coal event
            # NOTE: coal pos needs -.5 when using add_parent()

            tag, coal_time, coal_name, pos = event
            
            # find new_node lineage
            node1, parent1 = walk_up3(leaf, pos, coal_time, event_name)
            
            # find node that new_node coals with
            node2 = arg[coal_name]
            node2, parent2 = walk_up3(node2, pos, coal_time, event_name)

            if node1 == node2:
                print "SKIP"
                continue
            if parent1 and parent2:
                print "SKIP2"
                continue


            # create coal node
            coal = arg.new_node(age=coal_time, event="coal",
                                children=[node1, node2])
            new_nodes[coal] = event_name

            assert max(node1.age, node2.age) <= coal_time
            
            
            if parent1:
                # break branch
                node1.parents[node1.parents.index(parent1)] = coal
                parent1.children[parent1.children.index(node1)] = coal
                coal.parents.append(parent1)
                assert coal.age <= parent1.age
            else:
                # coal is off ARG
                add_parent(node1, coal, pos-.5)

            # is there a chance of the parent changing because of break with
            # node1 and parent1
            if parent2:
                # break branch
                node2.parents[node2.parents.index(parent2)] = coal
                parent2.children[parent2.children.index(node2)] = coal
                coal.parents.append(parent2)
                assert coal.age <= parent2.age, (coal.age, parent2.age)
            else:
                # coal is off ARG
                add_parent(node2, coal, pos-.5)

            nstubs -=1 

    for node in arg:
        for i, parent in enumerate(node.parents):
            if parent is None:
                print "NULL", node, i, node.event, node.pos, node.age

    


#=============================================================================
'''



'''

def get_deterministic_transitions2(states1, states2, times,
                                  tree, last_tree,
                                  recomb_branch, recomb_time,
                                  coal_branch, coal_time):

    # recomb_branch in tree
    # coal_branch in last_tree
    
    def find_state(node, time):
        b = util.INF
        state2 = None
        
        while len(node.children) == 1:
            node = node.children[0]
        
        for j, (n, t) in enumerate(states2):
            if node.name == n and time <= t < b:
                b = t
                state2 = j
        assert b == time
        assert state2 is not None, ((node, time), states2)
        return state2
    

    # get leaves under recomb_node
    recomb_leaves = set(last_tree.leaves(last_tree[recomb_branch]))
    
    next_states = []
    for i, state1 in enumerate(states1):
        node1, a = state1
        time = times[a]
        
        leaves1 = set(last_tree.leaves(last_tree[node1]))        
        remain = leaves1 - recomb_leaves
        

        if (node1, a) == (coal_branch, coal_time):
            # not a deterministic case (just mark i-->i)
            next_states.append(i)
        
        elif len(remain) > 0:
            # SPR only removes a subset of descendents
            # trace up from remaining leaf to find correct new state
            node = arg_lca(tree, [x.name for x in remain],
                           time, 0, ignore=recomb_branch)
            next_states.append(find_state(node, a))

        else:
            # SPR is on same branch as new chromosome
            if recomb_time >= a:
                # we move with SPR subtree
                node = arg_lca(tree, [x.name for x in leaves1], None, 0)
                next_states.append(find_state(node, a))

            else:
                # SPR should not be able to coal back onto same branch
                # this would be a self cycle
                if coal_branch == node1:
                    print (recomb_branch, recomb_time), \
                          (coal_branch, coal_time)
                    treelib.draw_tree_names(last_tree.get_tree(),
                                            minlen=8, maxlen=8)
                    treelib.draw_tree_names(tree.get_tree(),
                                            minlen=8, maxlen=8)

                    print "path1"
                    ptr = last_tree[recomb_branch]
                    ptr = ptr.parents[0]
                    while len(ptr.children) == 1:
                        print ptr.name, ptr.event
                        ptr = ptr.parents[0]

                    print "path2"
                    ptr = tree[recomb_branch]
                    ptr = ptr.parents[0]
                    while len(ptr.children) == 1:
                        print ptr.name, ptr.event
                        ptr = ptr.parents[0]
                    
                    assert False

                
                # SPR subtree moves out from underneath us
                # therefore therefore the new chromosome coalesces with
                # the branch above the subtree

                # search up for parent
                ptr = last_tree[recomb_branch]
                ptr = ptr.parents[0]
                while len(ptr.children) == 1:
                    ptr = ptr.parents[0]
                b = times.index(ptr.age)

                if ptr.name not in tree:
                    # we are above root
                    assert ptr.age >= tree.root.age
                    next_states.append(find_state(tree.root, b))
                else:
                    ptr = tree[ptr.name]
                    next_states.append(find_state(ptr, b))

    return next_states
'''

'''
def tree_lca(tree, leaves, time, pos, ignore=None, order=None):

    if order is None:
        order = dict((node, i) for i, node in enumerate(
            tree.postorder_marginal_tree(pos-.5)))
    local = set(order.keys())
    if ignore is not None and ignore in tree:
        ptr = tree[ignore]
        local.remove(ptr)
        ptr = ptr.parents[0] if ptr.parents else None

        while ptr and ptr in local:
            if len(ptr.children) == 2:
                break
            local.remove(ptr)
            #print "remove", ptr
            ptr = ptr.parents[0] if ptr.parents else None

    queue = [(order[tree[x]], tree[x]) for x in leaves]
    seen = set(x[1] for x in queue)
    heapq.heapify(queue)

    while len(queue) > 1:
        i, node = heapq.heappop(queue)
        parent = node.parents[0] if node.parents else None
        if parent and parent not in seen:
            seen.add(parent)
            heapq.heappush(queue, (order[parent], parent))
    node = queue[0][1]
    parent = node.parents[0] if node.parents else None

    # walk up appropriate time if given
    if time is not None:
        while parent and parent.age <= time:
            #print "up", parent
            if (len(parent.children) == 2 and
                parent.children[0] in local and
                parent.children[1] in local):
                #print "halt"
                break
            node = parent
            parent = node.parents[0] if node.parents else None

        if parent:
            if parent.age < time:
                print (leaves, parent.age, time, ignore)
                print local
                tree = tree.get_tree()
                treelib.draw_tree_names(tree, maxlen=8, minlen=8)
                treelib.remove_single_children(tree)
                tree.write()
                assert False

    return node
'''

