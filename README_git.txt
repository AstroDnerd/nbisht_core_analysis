
A procedure for git workflow purposes.

ON THE WEB
1. Get a Github account
2. Go to the project of interest and click on "Fork" on the top right hand
corner

ON YOUR TERMINAL
3. git clone URL_OF_FORK  # clone the fork
4. cd forkedRepo
5. git remote -v          # check current remotes, should see
"origin" and url of the fork if not do:
git remote add origin URL_OF_FORK

6. git remote add upstream URL_OF_PROJECT 
7. Repeat step 5 to check for the respective 'origin' and 'upstream'

# YOU ARE DONE SETTING UP!  ...now
# TO CONTRIBUTE, IFFFFF your group is small enough/advisor is nice enough
# 7-9 are best practices I've found without using a BRANCH, step 16

7. if doing this for the first time, be mindful of core .py documentation you
may have changed for your purpose of running things as the following will make
changes for everyone in your group! this is why maybe it's best to do step 11
first if multiple people are contributing...

8. make a directory with all of your .py files/work so as to separate your
contribution and to keep yourself organized
9. for steps 9 and 10 mv all of your files to that new directory, you could cp
these back when you're done updating the repos

10. git status            # it's nice to see what's going on
11. git add .             # if you want to add all the changes, otherwise follow
respecive prompts
12 git commit             # then it will prompt you to write a message
13. git pull upstream     # to sync with the project 
14. git push upstream     # to sync your contributions to the project
15. git push origin       # to sync your fork of the project

# OTHERWISE DO...
16. Go to:  https://www.dataschool.io/how-to-contribute-on-github/ 
    AND DO STEPS 8-18 to learn how to BRANCH

