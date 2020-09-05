library(devtools) 
install_github("HankStevens/primer", force=TRUE)
library(primer)


# to commit everything
# > git add -A && git commit -m 'staging all files'
# > git push -u origin master

# git add -A to add all files new files, changes and removed files.
# git commit -m  "Your message" to save the changes done in the files.
# git push -u origin master to send your committed changes to a remote repository, where the local branch is named master to the remote named origin

# /usr/bin/git push origin HEAD:refs/heads/master