# Git Workflow Guide for MATLAB Project

**Author:** JXLI
**Date:** 2025-02-09
**Purpose:** Version control for condensed matter physics tight-binding calculations

---

## 📚 Table of Contents

1. [Git Basics](#git-basics)
2. [Setting Up Git](#setting-up-git)
3. [Creating and Managing Branches](#creating-and-managing-branches)
4. [Making Commits](#making-commits)
5. [Working with Remote Repositories](#working-with-remote-repositories)
6. [Best Practices](#best-practices)
7. [Common Workflows](#common-workflows)
8. [Troubleshooting](#troubleshooting)

---

## 🎯 Git Basics

### What is Git?

Git is a **distributed version control system** that tracks changes in your code over time. It allows you to:
- Save snapshots of your work (commits)
- Create parallel development lines (branches)
- Collaborate with others
- Revert to previous versions
- Track who changed what and when

### Key Concepts

- **Repository (Repo):** A directory tracked by Git
- **Commit:** A snapshot of your code at a point in time
- **Branch:** An independent line of development
- **Remote:** A version of your repository hosted elsewhere (e.g., GitHub, GitLab)
- **Clone:** Download a copy of a remote repository
- **Pull:** Download changes from remote
- **Push:** Upload your changes to remote

---

## 🚀 Setting Up Git

### 1. Check Git Installation

```bash
git --version
```

If not installed, download from [git-scm.com](https://git-scm.com/)

### 2. Configure Git (First Time Only)

```bash
# Set your name
git config --global user.name "Your Name"

# Set your email
git config --global user.email "your.email@example.com"

# Set default editor (optional)
git config --global core.editor "vim"

# View your configuration
git config --list
```

### 3. Initialize a New Repository

```bash
# Navigate to your project directory
cd /Volumes/T9/work/tb/matlab/matlab_organized_v1/

# Initialize Git repository
git init

# Check status
git status
```

**Note:** The parent directory is already a Git repository, so we'll work within that.

---

## 🌿 Creating and Managing Branches

### Understanding Branches

Branches allow you to work on different features or experiments without affecting the main code. Think of it like parallel universes of your code.

```
main ──────●──────●──────●──────●
            \
             \── feature ──●──●──●
```

### Current Repository Status

```bash
# Check current branch
git branch

# Output: * main  (you're on main branch)
```

### Creating a New Branch

```bash
# Create and switch to new branch (recommended)
git checkout -b organize-v1

# Or in two steps:
git branch organize-v1      # Create branch
git checkout organize-v1    # Switch to branch

# Verify you're on the new branch
git branch
# Output:
#   main
# * organize-v1
```

### Branch Naming Conventions

Good branch names:
- `feature/quantum-geometry` - New feature
- `fix/berry-curvature-bug` - Bug fix
- `organize-v1` - Project reorganization
- `experiment/new-algorithm` - Experimental work

Bad branch names:
- `test` - Too vague
- `branch1` - Not descriptive
- `my-branch` - Not informative

### Switching Between Branches

```bash
# Switch to main branch
git checkout main

# Switch back to organize-v1
git checkout organize-v1

# Create and switch in one command (modern Git)
git switch -c organize-v1
```

### Listing Branches

```bash
# List local branches
git branch

# List all branches (including remote)
git branch -a

# List with last commit info
git branch -v
```

### Deleting Branches

```bash
# Delete a branch (safe - won't delete if unmerged)
git branch -d branch-name

# Force delete (careful!)
git branch -D branch-name
```

---

## 💾 Making Commits

### The Git Workflow

```
Working Directory → Staging Area → Repository
     (edit)           (git add)     (git commit)
```

### 1. Create .gitignore First

```bash
# Navigate to project root
cd /Volumes/T9/work/tb/matlab/matlab_organized_v1/

# Create .gitignore file
cat > .gitignore << 'EOF'
# MATLAB temporary files
*.asv
*.m~
*.swp

# MATLAB data files (usually too large)
*.mat
*.dat
*.h5
*.hdf5

# Images and figures
*.png
*.jpg
*.jpeg
*.eps
*.pdf
*.fig
*.tiff

# Compressed files
*.zip
*.tgz
*.tar.gz
*.rar

# System files
.DS_Store
Thumbs.db
*.tmp

# Data directory (symlink or large data)
data/
*/data/

# Build artifacts
*.mex*
*.o
*.a

# IDE files
.vscode/
.idea/
*.sublime-*

# Backup files
*~
*.bak
*.backup

# Log files
*.log

# Keep structure but ignore specific large files
# You can add specific exceptions with !pattern
EOF
```

### 2. Check What Will Be Committed

```bash
# See status
git status

# See detailed changes
git diff

# See which files are new
git status -u
```

### 3. Stage Files

```bash
# Stage all files
git add .

# Stage specific file
git add filename.m

# Stage multiple files
git add file1.m file2.m

# Stage all .m files
git add "*.m"

# Stage by directory
git add +MTB/

# Interactive staging (choose what to stage)
git add -i
```

### 4. Review Staged Changes

```bash
# See what's staged
git diff --staged

# Or
git diff --cached
```

### 5. Make a Commit

```bash
# Commit with message
git commit -m "Initial commit: organized project structure"

# Commit with detailed message
git commit -m "Add quantum geometry module

- Implemented quantum metric calculation
- Added Berry curvature dipole functions
- Updated documentation"

# Commit all tracked changes (skip staging)
git commit -am "Quick fix: typo in README"
```

### Good Commit Messages

**Format:**
```
<type>: <subject>

<body>

<footer>
```

**Types:**
- `feat:` New feature
- `fix:` Bug fix
- `docs:` Documentation
- `refactor:` Code restructuring
- `test:` Adding tests
- `chore:` Maintenance

**Examples:**

```bash
git commit -m "feat: add TaIrTe4 quantum metric module"

git commit -m "fix: correct Berry curvature calculation in get_bcd.m"

git commit -m "docs: update MTB toolbox README with examples"

git commit -m "refactor: reorganize project by material systems"
```

### Viewing Commit History

```bash
# View commit log
git log

# One line per commit
git log --oneline

# Show changes in each commit
git log -p

# Last N commits
git log -n 5

# Graphical view
git log --graph --oneline --all

# Search commits
git log --grep="quantum"
```

---

## 🌐 Working with Remote Repositories

### Understanding Remotes

A **remote** is a version of your repository hosted on a server (GitHub, GitLab, etc.).

### Popular Git Hosting Services

1. **GitHub** (github.com) - Most popular
2. **GitLab** (gitlab.com) - Good for private repos
3. **Bitbucket** (bitbucket.org) - Atlassian product

### Setting Up Remote

#### Option 1: Clone Existing Repository

```bash
# Clone from GitHub
git clone https://github.com/username/repository.git

# Clone with specific branch
git clone -b branch-name https://github.com/username/repository.git
```

#### Option 2: Add Remote to Existing Repo

```bash
# Add remote
git remote add origin https://github.com/username/matlab-tb.git

# Verify remote
git remote -v

# Output:
# origin  https://github.com/username/matlab-tb.git (fetch)
# origin  https://github.com/username/matlab-tb.git (push)
```

### Pushing to Remote

```bash
# Push current branch to remote
git push origin organize-v1

# Push and set upstream (first time)
git push -u origin organize-v1

# After setting upstream, just:
git push

# Push all branches
git push --all origin

# Force push (DANGEROUS - use carefully!)
git push --force origin organize-v1
```

### Pulling from Remote

```bash
# Fetch and merge changes
git pull origin main

# Just fetch (don't merge)
git fetch origin

# Pull specific branch
git pull origin organize-v1
```

### Creating a GitHub Repository

1. Go to [github.com](https://github.com) and sign in
2. Click "+" → "New repository"
3. Name: `matlab-condensed-matter`
4. Description: "Tight-binding calculations for condensed matter physics"
5. Choose Public or Private
6. Don't initialize with README (we have one)
7. Click "Create repository"
8. Follow instructions to push existing repo:

```bash
git remote add origin https://github.com/yourusername/matlab-condensed-matter.git
git branch -M main
git push -u origin main
```

---

## ✅ Best Practices

### 1. Commit Frequently

- Commit logical units of work
- Don't wait until end of day
- Small commits are easier to review and revert

### 2. Write Meaningful Commit Messages

❌ Bad:
```bash
git commit -m "update"
git commit -m "fix stuff"
git commit -m "asdfasdf"
```

✅ Good:
```bash
git commit -m "feat: implement Wilson loop for TaIrTe4"
git commit -m "fix: correct k-point mesh generation in HFMF"
git commit -m "docs: add usage examples for quantum geometry"
```

### 3. Use Branches

- `main` or `master`: Stable production code
- `develop`: Active development
- `feature/*`: New features
- `fix/*`: Bug fixes
- `experiment/*`: Experimental work

### 4. Keep .gitignore Updated

Don't commit:
- Large data files (*.mat, *.dat)
- Generated files (*.png, *.pdf)
- Temporary files (*.asv, *.m~)
- System files (.DS_Store)

### 5. Pull Before Push

```bash
git pull origin main
git push origin main
```

### 6. Review Changes Before Committing

```bash
git status
git diff
git diff --staged
```

---

## 🔄 Common Workflows

### Workflow 1: Feature Development

```bash
# 1. Create feature branch
git checkout -b feature/new-algorithm

# 2. Make changes
# ... edit files ...

# 3. Stage and commit
git add .
git commit -m "feat: implement new HFMF algorithm"

# 4. Push to remote
git push -u origin feature/new-algorithm

# 5. Create Pull Request on GitHub
# ... review and merge ...

# 6. Switch back to main and update
git checkout main
git pull origin main

# 7. Delete feature branch
git branch -d feature/new-algorithm
```

### Workflow 2: Bug Fix

```bash
# 1. Create fix branch
git checkout -b fix/berry-curvature

# 2. Fix the bug
# ... edit files ...

# 3. Commit
git commit -am "fix: correct sign error in Berry curvature"

# 4. Push and merge
git push origin fix/berry-curvature
```

### Workflow 3: Experiment

```bash
# 1. Create experiment branch
git checkout -b experiment/fft-acceleration

# 2. Try new approach
# ... experimental code ...

# 3. If successful, commit
git commit -am "experiment: FFT acceleration works!"

# 4. Merge to main or develop
git checkout develop
git merge experiment/fft-acceleration

# 5. If failed, just delete branch
git branch -D experiment/fft-acceleration
```

### Workflow 4: Daily Work

```bash
# Morning: Update your code
git checkout main
git pull origin main

# Create daily work branch
git checkout -b work/2025-02-09

# Work on your code
# ... edit, add, commit multiple times ...

# End of day: Push your work
git push origin work/2025-02-09
```

---

## 🆘 Troubleshooting

### Problem 1: Accidentally Committed Wrong Files

```bash
# Remove from staging (before commit)
git reset HEAD filename.m

# Undo last commit (keep changes)
git reset --soft HEAD~1

# Undo last commit (discard changes) - CAREFUL!
git reset --hard HEAD~1

# Amend last commit
git commit --amend -m "corrected commit message"
```

### Problem 2: Want to Discard Local Changes

```bash
# Discard changes to specific file
git checkout -- filename.m

# Discard all changes
git reset --hard HEAD

# Clean untracked files
git clean -fd
```

### Problem 3: Merge Conflicts

```bash
# During merge, if conflicts occur:
git merge feature-branch

# Output: CONFLICT (content): Merge conflict in file.m

# 1. Open file.m and look for:
<<<<<<< HEAD
your code
=======
their code
>>>>>>> feature-branch

# 2. Edit to resolve conflict

# 3. Stage resolved file
git add file.m

# 4. Complete merge
git commit -m "merge: resolve conflicts in file.m"
```

### Problem 4: Need to Stash Changes

```bash
# Save work in progress
git stash

# Do other work (switch branches, pull, etc.)
git checkout main
git pull

# Restore your work
git stash pop

# List stashes
git stash list

# Apply specific stash
git stash apply stash@{0}
```

### Problem 5: Large Files Committed by Mistake

```bash
# Remove from Git but keep locally
git rm --cached large_file.mat

# Remove from history (advanced)
git filter-branch --force --index-filter \
  'git rm --cached --ignore-unmatch large_file.mat' \
  --prune-empty --tag-name-filter cat -- --all
```

---

## 📋 Quick Reference Cheat Sheet

### Basic Commands

```bash
git status              # Check status
git add .               # Stage all changes
git commit -m "msg"     # Commit with message
git push                # Push to remote
git pull                # Pull from remote
```

### Branch Commands

```bash
git branch              # List branches
git checkout -b name    # Create and switch to branch
git checkout name       # Switch to branch
git merge name          # Merge branch into current
git branch -d name      # Delete branch
```

### History Commands

```bash
git log                 # View history
git log --oneline       # Compact history
git diff                # Show changes
git show commit-id      # Show specific commit
```

### Remote Commands

```bash
git remote -v           # List remotes
git remote add origin URL  # Add remote
git push origin branch  # Push branch
git pull origin branch  # Pull branch
git clone URL           # Clone repository
```

### Undo Commands

```bash
git reset HEAD file     # Unstage file
git checkout -- file    # Discard changes
git revert commit-id    # Revert commit
git reset --hard HEAD~1 # Undo last commit (DANGER!)
```

---

## 🎓 Learning Resources

### Tutorials
- [Official Git Tutorial](https://git-scm.com/docs/gittutorial)
- [GitHub Guides](https://guides.github.com/)
- [Atlassian Git Tutorial](https://www.atlassian.com/git/tutorials)

### Interactive Learning
- [Learn Git Branching](https://learngitbranching.js.org/)
- [Git-it Tutorial](https://github.com/jlord/git-it-electron)

### Books
- [Pro Git](https://git-scm.com/book/en/v2) (Free online)
- [Git Pocket Guide](https://www.oreilly.com/library/view/git-pocket-guide/9781449327507/)

### Cheat Sheets
- [GitHub Git Cheat Sheet](https://education.github.com/git-cheat-sheet-education.pdf)
- [Atlassian Git Cheat Sheet](https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet)

---

## 🚀 Next Steps for This Project

### Immediate Actions

1. **Create .gitignore**
   ```bash
   cd /Volumes/T9/work/tb/matlab/matlab_organized_v1/
   # Use the .gitignore template provided above
   ```

2. **Create Branch**
   ```bash
   git checkout -b organize-v1
   ```

3. **Initial Commit**
   ```bash
   git add .
   git commit -m "Initial commit: organized MATLAB TB project

   - Created organized project structure
   - Added MTB toolbox with 90+ functions
   - Added tbHFMF module
   - Included all material system projects
   - Created comprehensive documentation"
   ```

4. **Push to Remote** (after setting up GitHub)
   ```bash
   git push -u origin organize-v1
   ```

### Future Workflow

```bash
# Daily workflow
git checkout organize-v1
git pull origin organize-v1

# Make changes
# ... work on code ...

git add .
git commit -m "feat: add TaIrTe4 quantum geometry examples"
git push

# When ready to merge to main
git checkout main
git merge organize-v1
git push origin main
```

---

**Last Updated:** 2025-02-09
**Version:** 1.0
**Maintained by:** JXLI
