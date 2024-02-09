
A procedure for git workflow purposes.

ON THE WEB
1. Get a Github account
2. Go to the project of interest and click on "Fork" on the top right hand
corner

ON YOUR TERMINAL
3. git clone URL_OF_FORK  # clone the fork, USING SSH!! (https=pain)
4. cd forkedRepo
5. git remote -v          # check current remotes, should see
"origin" and url of the fork if not do:
git remote add origin URL_OF_FORK

6. git remote add upstream URL_OF_PROJECT 
7. Repeat step 5 to check for the respective 'origin' and 'upstream'

# YOU ARE DONE SETTING UP!  ...now
# TO CONTRIBUTE, IFFFFF your group is small enough/advisor is nice enough
# 7-9 are best practices I've found without using a BRANCH

8. if doing this for the first time, be mindful of core .py documentation you
may have changed for your purpose of running things as the following will make
changes for everyone in your group! this is why maybe it's best to do step 16
first if multiple people are contributing...

9. make a directory with all of your .py files/work so as to separate your
contribution and to keep yourself organized

10. git status            # it's nice to see what you have changed 
11. git add .             # if you want to add all the changes, otherwise only
add the files you want

12 git commit -m 'message of change'   # a brief message of changes 
13. git pull upstream     # to sync with the project BEFORE pushing

14. git push upstream     # to sync your contributions to the fork of your
project, may need an extra word after upstream
15. git push origin       # to sync your fork to the original project, may need
an extra word after origin

# OTHERWISE DO...
16. Go to:  https://www.dataschool.io/how-to-contribute-on-github/ 
    AND DO STEPS 8-18 to learn how to BRANCH

# OR the simples procedure after #7:
8. git checkout -b 'branchname'  # however, it may be enough to have forked 
9. Repeat steps 10-15, except that it is not necessary to add anything post pull/push.
10. From github (web), pick the branch & create a pull request (PR)
