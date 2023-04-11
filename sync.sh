

# go to sync branch
git checkout sync

# to update to the dev branch in the current repo
git fetch
git rebase dev


# to get latest updates from GitHub
git checkout sync
git remote add upstream https://github.com/TRON-Bioinformatics/splice2neo.git
git status


# resolve conflicts manually, e.g. fix version in DESCRIPTION
git push

