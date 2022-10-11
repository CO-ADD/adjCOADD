# Get start with Git
This document includes two parts: using Git as code version control locally and collabrating with remote project




## Collabrating with remote project on github.						
1. using ```git clone the repo url``` for the first time. And after that every time using ```git pull origin main``` or ```git pull origin master``` before doing any job.						
						
2. before doing you job leave the main/master branch check into another branch using ```git checkout -b <branchname>```						
						
3. after your work, still in the same branch push to remote github with ```git push origin <branchname>```					
						
4. goto remote project on github. Your changing record will show on the top and next to it is "compare and pull request" green button. Click the button.						
						
5. using "Review" in right bar to request a review and approve.						
						
Note: NOT do any job directly on the local main/master branch. And it should be kept updated before your start to work. 						
						
refer to : https://www.youtube.com/watch?v=MnUd31TvBoU						


## locally using git

For a collabration you can skip 1.-3. and 10.

1. goto your django project root directory: 
```
$: cd myDjangoproject
/myDjangoproject $:
```
2. initialize git repo with git init:
```
/myproject $: git init -b main (or master)
```
this will create a .git/ directory in your project root to save your code in different version. And you will see the information :" Initialized empty Git repository in your/projcet path/.git/..." or similar…

3. create .gitignore file in your project root. You can copy anyone's in github or somewhere… (for python or django).This file can be used to list the files that you do not want git to track

4. using ```git add .``` or ```git add <filesname> ```add changed codes files to the staging area for Git

5. ```git commit``` or ```git commit -m "message" ```  record the changes made to the files to a local repository. For easy reference, each commit has a unique ID(you can use it to roll back to this version).

6. make changes:Unstage Changes for a <file> (After git add, Before git commit)
use ```git reset HEAD <file>...```to unstage
use ```git checkout -- <file>...``` to discard changes in working directory

7. Undo Commit (After git commit, Before git push)
the same like 6S

8. check repo state using ```git status```
10. create an empty repo in github, add the remote url to git using ```git remote add origin <REMOTE_URL>``` then push your project to it with ``` git push origin main```.