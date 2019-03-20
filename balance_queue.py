"""
### note
# right now this isn't using the RAC, and will exit if there are not more than one def-X accounts
###

###
# usage: balance_queue.py phaseOFpipeline
###

###
# purpose: evenly redistributes jobs across available accounts based on priority based on the job name
#          (phaseOFpipeline);
#          helps speed up effective run time
###
"""

import os, shutil, sys, math, subprocess


def announceacctlens(accounts, fin):
    print('%s job announcement' % ('final' if fin is True else 'first'))
    for account in accounts:
        print('%s jobs on %s' % (str(len(accounts[account])), account))


def checksq(sq):
    exitneeded = False
    if not isinstance(sq, list):
        print("type(sq) != list, exiting %(thisfile)s" % globals())
        exitneeded = True
#     if len(sq) == 0:  # i don't think I need this with the new -h in subproc in getsq()
#         print("len(sq) == 0, exiting %(thisfile)s" % globals())
#         exitneeded = True
    for s in sq:
        if not s == '':
            if 'socket' in s.lower():
                print("socket in sq return, exiting %(thisfile)s" % globals())
                exitneeded = True
            if not int(s.split()[0]) == float(s.split()[0]):
                print("could not assert int == float, %s" % (s[0]))
                exitneeded = True
    if exitneeded is True:
        print('slurm screwed something up for %(thisfile)s, lame' % globals())
        exit()
    else:
        return sq


def getsq(grepping, states=[], balance=False):
    if isinstance(grepping, str):
        # in case I pass a single str instead of a list of strings
        grepping = [grepping]

    # get the queue, without a header
    sq = subprocess.check_output([shutil.which('squeue'),
                                  '-u',
                                  os.environ['USER'],
                                  '-h']).decode('utf-8').split('\n')
    sq = [s for s in sq if not s == '']
    checksq(sq)  # make sure slurm gave me something useful

    # look for the things I want to grep (serial subprocess.Popen() are a pain with grep)
    grepped = []
    if len(sq) > 0:
        for q in sq:  # for each job in queue
            splits = q.split()
            if 'CG' not in splits:  # grep -v 'CG'
                keepit = 0
                for split in splits:
                    for grep in grepping:  # see if all necessary greps are in the job
                        if grep.lower() in split.lower():
                            keepit += 1
                if keepit == len(grepping):
                    # see if any of the state conditions are met (does not need all states, obviously)
                    if len(states) > 0:
                        print(states, splits[4])
                        keepit2 = False
                        for state in states:
                            if state.lower() == splits[4].lower():
                                keepit2 = True
                        if keepit2 is True:
                            grepped.append(tuple(splits))
                    else:
                        grepped.append(tuple(splits))

        if len(grepped) > 0:
            return grepped
    else:
        print('no jobs in queue')
        if balance is True:
            exit()
        else:
            return sq


def adjustjob(acct, jobid):
    print(acct)
    subprocess.Popen([shutil.which('scontrol'), 'update', 'Account=%s' % acct, 'JobId=%s' % str(jobid)])
    # os.system('scontrol update Account=%s_cpu JobId=%s' % (acct, str(jobid)) )


def getaccounts(sq, stage):
    accounts = {}
    for q in sq:
        pid = q[0]
        account = q[2]
        account = account.split("_")[0]
        if account not in accounts:
            accounts[account] = {}
        accounts[account][pid] = q
#     if len(accounts.keys()) == 3 and stage != 'final': # all accounts have low priority ### use 3 when using RAC
    if len(accounts.keys()) == 2 and stage != 'final':  # all accounts have low priority   ### use 2 when not using RAC
        print('all accounts have low priority, leaving queue as-is')
        announceacctlens(accounts, True)
        exit()
    return accounts


def getbalance(accounts, num):
    sums = 0
    for account in accounts:
        sums += len(accounts[account].keys())
    bal = math.ceil(sums/num)
    print('bal%i %i= ' % (num, bal))
    return bal


# def checknumaccts(accts, checking, mc):
#     # len(accounts) will never == 2 after pop, since I checked for len(accounts) == 3
#     if len(accts.keys()) == 0:
#         if checking == 'RAC':
#             print('RAC has low priority status, skipping RAC as taker')
#         else:
#             print('moved %s jobs to RAC' % str(mc))
#         exit()
#
#
# def redistribute4g(accounts, bal, rac, mcount=0):
#     if rac in accounts:    # no need to redistribute to rac if rac has low priority
#         accounts.pop(rac)  # drop rac from list to redistribute, exit if nothing to redistribute
#         checknumaccts(accounts, 'rac', '')    # if all jobs are on rac, exit
#         return accounts
#     keys = list(accounts.keys())
#     print('before loop %s' % keys)
#     for account in keys:
#         # distribute 4G jobs to rac
#         pids = list(accounts[account].keys())
#         mcount = 0
#         for pid in pids:
#             mem = int([m for m in accounts[account][pid] if m.endswith('M')][0].split("M")[0])
#             if mem <= 4000:
#                 # if it can be scheduled on the rac, change the account of the jobid, and remove jobid from list
#                 adjustjob(rac, pid)
#                 accounts[account].pop(pid)
#                 mcount += 1
#                 if mcount == bal:
#                     break
#         print("distributed {} jobs from {} to rac".format(mcount, account))
#         if len(accounts[account].keys()) == 0:
#             accounts.pop(account)
#     checknumaccts(accounts, 'none', mcount)  # if all jobs were redistributed to the rac, exit
#     return accounts


def gettaker(accounts, defs):
    print('gettaker defs =', defs)
    giver = ''
    keys = list(accounts.keys())
    if len(keys) == 2:
        # if there are two accounts, figure out which account has more
        maxx = 0
        for acct in keys:
            if len(accounts[acct]) > maxx:
                giver = acct
                maxx = len(accounts[acct])
    else:
        if not len(keys) == 1:
            print('assertion error')
        giver = keys[0]
    # taker = list('def-saitken', 'def-yeaman'}.symmetric_difference({giver}))[0]
    taker = list(set(defs).symmetric_difference({giver}))[0]
    return giver, taker


def givetotaker(giver, taker, accounts, bal):
    taken = 0
    pids = list(accounts[giver].keys())
    numtotake = len(pids) - bal
    if bal == 1 and len(pids) == 1:
        numtotake = 1
    printout = 'giver has {} jobs to give. (bal= {}). Giver ({}) is giving {} jobs to taker ({})'.format(len(pids),
                                                                                                         bal,
                                                                                                         giver,
                                                                                                         numtotake,
                                                                                                         taker)
    print("\t %s" % printout)
    if numtotake > 0:
        for pid in pids[::-1]:  # re-assign the newer jobs, hopefully older jobs will eventually run
            adjustjob(taker, pid)
            taken += 1
            if taken == numtotake:
                print("\t redistributed %s jobs from %s to %s" % (str(taken), giver, taker))
                break
    else:
        print("\t giver sees that taker has enough, so giver is not giving")


def get_availaccounts():
    accts = subprocess.check_output([shutil.which('sshare'),
                                     '-U',
                                     '--user',
                                     os.environ['USER'],
                                     '--format=Account']).decode('utf-8').split('\n')
    accts = [acct.split()[0] for acct in accts if 'def' and 'cpu' in acct]
    defs = [acct for acct in accts if 'def' in acct]
    rac = [acct for acct in accts if 'rrg' in acct]
    if len(defs) == 1:
        # no need to try and balance
        print(f'there is only one account ({defs[0]}), no need to balance queue.\nexiting balance_queue.py')
        exit()
    if len(rac) == 1:
        rac = rac[0]
    elif len(rac) == 0:
        rac = ''
    print('accts = ', accts)
    return defs, rac


def main(thisfile, phase):
    globals().update({'thisfile': thisfile, 'phase': phase})

    # get accounts available for billing
    defs, rac = get_availaccounts()

    # get the queue
    sq = getsq(grepping=[phase, 'Priority'], balance=True)

    # get per-account counts of jobs in Priority pending status, exit if all accounts have low priority
    accts = getaccounts(sq, '')
    announceacctlens(accts, False)

    # figure out how many to balance remaining
    # balance = getbalance(accts, 3)

    # redistribute 4G jobs to RAC unless RAC has low priority, exit if all jobs redistributed or no jobs to redistribute
    # accts = redistribute4g(accts, balance, rac)

    # figure out which account to add to
    giver, taker = gettaker(accts, defs)

    # redistribute to taker
    balance = getbalance(accts, 2)
    givetotaker(giver, taker, accts, balance)

    # announce final job counts
    announceacctlens(getaccounts(getsq(grepping=[phase, 'Priority']),
                                 'final'),
                     True)


if __name__ == '__main__':
    # args
    thisfile, phase = sys.argv

    main(thisfile, phase)
