# BIOSCIENCE PYTHON TUTORIALS
### For UW Bioengineering 2019-2020 honors

ERIC AND MAX:
```
git clone <repo_url>
git checkout -b <branch_name>
```
code, add and commit your changes...
```
git add -p
git commit -m "<your message>"
```
Continue adding and commiting to log and save your changes. When you're done with small changes:
```
git push origin <branch_name>
```
When you're done with your module:
* go to github --> branches --> create pull request.
* wait for everyone to review and approve
* merge

Now your changes are on master!

Clean up your local repo and pull new changes:
```
git checkout master
git pull
git branch -d <branch_name>
```

Then to start your next module:
```
git checkout -b <next_branch_name>
```
