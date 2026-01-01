# Git Worktrees Guide - Understanding Your Setup

**Last Updated:** 2025-01-14  
**Purpose:** Complete guide to Git worktrees, branch switching, and Cursor integration

---

## 🎯 **WHAT ARE GIT WORKTREES?**

**Git worktrees** allow you to have **multiple working directories** for the **same repository**, each checking out different branches or commits simultaneously.

### **Traditional Git (Single Working Directory)**
```
crispr-assistant-main/
├── .git/              # Repository metadata
├── src/               # Your code
└── ...                # You can only work on ONE branch at a time
```

**Problem:** To switch branches, you must:
1. Commit or stash your changes
2. `git checkout other-branch`
3. Work on that branch
4. Switch back (repeat)

### **Git Worktrees (Multiple Working Directories)**
```
Main Repo:
/Users/fahadkiani/Desktop/development/crispr-assistant-main/
├── .git/             # Main repository
├── src/                 # Working on 'main' branch
└── ...

Worktree 1:
/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw/
├── .git -> (points to main repo's .git/worktrees/nkw)
├── src/               # Working on 'latest-backend-updates' branch
└── ...

Worktree 2:
/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/rxk/
├── .git -> (points to main repo's .git/worktrees/rxk)
├── src/               # Working on another branch/commit
└── ...
```

**Benefit:** Work on **multiple branches simultaneously** without switching!

---

## 📊 **YOUR CURRENT SETUP**

### **Main Repository**
- **Location:** `/Users/fahadkiani/Desktop/development/crispr-assistant-main`
- **Current Branch:** `main`
- **Commit:** `e830b5c` (HEAD)
- **Remotes:**
  - `origin`: `git@github.com:fjkiani/crispro.git`
  - `upstream`: `https://github.com/ArcInstitute/evo2.git`

### **Worktree 1: `nkw`**
- **Location:** `/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw`
- **Status:** Detached HEAD at `e830b5c`
- **Purpose:** Cursor's working directory (this is where you're currently working)
- **Changes:** 
  - Modified: `.cursor/GITHUB_REPOSITORIES.md`
  - Modified: `.specstory/history/2026-01-01_05-56Z-git-repository-context.md`
  - Untracked: `delete_all_repos.sh`, `delete_repos.sh`

### **Worktree 2: `rxk`**
- **Location:** `/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/rxk`
- **Status:** Detached HEAD at `e830b5c`
- **Purpose:** Another Cursor worktree (possibly for a different session)

---

## 🔗 **HOW WORKTREES ARE LINKED**

### **The `.git` File in Worktrees**

In a worktree, the `.git` file is **not a directory** - it's a **text file** that points to the main repository:

```bash
# In worktree: /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw/.git
gitdir: /Users/fahadkiani/Desktop/development/crispr-assistant-main/.git/worktrees/nkw
```

This tells Git: "The actual `.git` directory is over there in the main repo's worktrees folder."

### **Shared Repository Data**

All worktrees share:
- ✅ **Same commit history** (all branches, all commits)
- ✅ **Same remotes** (origin, upstream)
- ✅ **Same refs** (branch pointers, tags)
- ✅ **Same config** (user.name, user.email, etc.)

Each worktree has:
- ✅ **Independent working directory** (different files checked out)
- ✅ **Independent index** (different staged changes)
- ✅ **Independent HEAD** (can be on different branches/commits)

---

## 🔄 **SWITCHING BRANCHES IN CURSOR**

### **Method 1: Switch Branch in Current Worktree**

**Current State:** Your worktree (`nkw`) is in **detached HEAD** state (not on a branch).

**To switch to a branch:**

```bash
# In Cursor terminal (or command palette)
cd /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw

# Option A: Create a new branch from current commit
git checkout -b my-new-branch

# Option B: Switch to existing branch (if it exists)
git checkout main
git checkout latest-backend-updates

# Option C: Create branch and track remote
git checkout -b feature-branch origin/main
```

**After switching:** Your worktree will be on that branch, and you can commit/push normally.

---

### **Method 2: Create a New Worktree for a Branch**

**Use this when:** You want to work on multiple branches **simultaneously** without switching.

```bash
# From main repo
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main

# Create new worktree for a branch
git worktree add ../worktrees/feature-branch feature-branch

# Or create worktree and new branch at same time
git worktree add ../worktrees/new-feature -b new-feature
```

**Result:** New folder created:** `/Users/fahadkiani/Desktop/development/worktrees/feature-branch/` with that branch checked out.

**In Cursor:** You can open this folder as a separate workspace!

---

### **Method 3: Use Cursor's Git UI**

**Cursor's Source Control Panel:**
1. Click the **Source Control** icon in sidebar (or `Cmd+Shift+G`)
2. Click the **branch name** at bottom-left (shows current branch/HEAD)
3. Select a branch from the list, or type to create new branch

**Cursor's Command Palette:**
1. Press `Cmd+Shift+P` (Mac) or `Ctrl+Shift+P` (Windows/Linux)
2. Type: `Git: Checkout to...`
3. Select branch or create new one

---

## 🎯 **BEST PRACTICES FOR WORKTREES**

### **1. Main Repo = Production/Stable Branch**

**Recommendation:** Keep `main` branch in the main repo directory.

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
git checkout main
# Work on stable/production code here
```

### **2. Worktrees = Feature Branches**

**Recommendation:** Use worktrees for feature development.

```bash
# Create worktree for feature
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
git worktree add ../worktrees/feature-xyz feature-xyz

# Open in Cursor
# File → Open Folder → /Users/fahadkiani/Desktop/development/worktrees/feature-xyz
```

### **3. Cursor Worktrees = Temporary/Experimental**

**Current Setup:** Cursor creates worktrees automatically (like `nkw`, `rxk`).

**Recommendation:** 
- Use Cursor worktrees for **quick experiments**
- When done, **merge back to main** and **remove worktree**

```bash
# After merging feature
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
git worktree remove ../worktrees/feature-xyz
```

---

## 🛠️ **COMMON WORKTREE COMMANDS**

### **List All Worktrees**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
git worktree list
```

**Output:**
```
/Users/fahadkiani/Desktop/development/crispr-assistant-main    e830b5c [main]
/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw  e830b5c (detached HEAD)
/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/rxk  e830b5c (detached HEAD)
```

### **Create New Worktree**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
git worktree add ../worktrees/my-feature my-feature-branch
```

### **Remove Worktree**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
git worktree remove ../worktrees/my-feature
# Or use --force if worktree has uncommitted changes
git worktree remove --force ../worktrees/my-feature
```

### **Move Worktree**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
git worktree move ../worktrees/old-name ../worktrees/new-name
```

### **Prune Orphaned Worktrees**
```bash
# Remove worktrees that no longer exist on disk
git worktree prune
```

---

## ⚠️ **DETACHED HEAD STATE - WHAT IT MEANS**

### **What is Detached HEAD?**

**Detached HEAD** means you're not on a branch - you're on a specific commit.

**Your Current State:**
```bash
cd /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw
git status
# Output: "Not currently on any branch."
```

### **Why Cursor Creates Detached HEAD Worktrees**

Cursor may create worktrees in detached HEAD state to:
- Avoid branch conflicts
- Work on specific commits
- Isolate experimental changes

### **How to Fix Detached HEAD**

**Option 1: Create Branch from Current Commit**
```bash
cd /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw
git checkout -b my-branch-name
# Now you're on a branch!
```

**Option 2: Switch to Existing Branch**
```bash
cd /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw
git checkout main
# Now you're on main branch
```

**Option 3: Create Branch and Push**
```bash
cd /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw
git checkout -b feature-branch
git push -u origin feature-branch
```

---

## 🔍 **UNDERSTANDING YOUR CURRENT SITUATION**

### **Current Worktree Status**

**Location:** `/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw`

**Status:**
- ✅ **Detached HEAD** at commit `e830b5c`
- ✅ **Has uncommitted changes:**
  - Modified: `.cursor/GITHUB_REPOSITORIES.md`
  - Modified: `.specstory/history/2026-01-01_05-56Z-git-repository-context.md`
  - Untracked: `delete_all_repos.sh`, `delete_repos.sh`

### **What to Do Next**

**Option A: Commit Changes to New Branch**
```bash
cd /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw
git checkout -b github-repo-cleanup
git add .
git commit -m "Add GitHub repository organization and cleanup scripts"
git push -u origin github-repo-cleanup
```

**Option B: Commit to Main Branch**
```bash
cd /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw
git checkout main
git add .
git commit -m "Add GitHub repository organization and cleanup scripts"
git push origin main
```

**Option C: Stash Changes and Switch**
```bash
cd /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw
git stash
git checkout main
# Work on main, then later:
git stash pop
```

---

## 🎨 **CURSOR-SPECIFIC WORKFLOW**

### **Opening Different Worktrees in Cursor**

**Method 1: File → Open Folder**
1. `File` → `Open Folder...`
2. Navigate to worktree directory
3. Cursor opens it as a new workspace

**Method 2: Command Palette**
1. `Cmd+Shift+P` → `File: Open Folder...`
2. Select worktree directory

**Method 3: Multiple Windows**
- Each worktree can be opened in a **separate Cursor window**
- Work on multiple branches simultaneously!

### **Cursor's Git Integration**

**Source Control Panel:**
- Shows changes in **current worktree**
- Commits affect **current worktree's branch/HEAD**
- Push/pull affects **current worktree's remote**

**Branch Switcher:**
- Bottom-left shows current branch/HEAD
- Click to switch branches (in current worktree)

---

## 📋 **QUICK REFERENCE**

### **Check Current Worktree**
```bash
git worktree list
```

### **Check Current Branch/HEAD**
```bash
git branch --show-current
# Or if detached HEAD:
git rev-parse --abbrev-ref HEAD
```

### **See All Branches**
```bash
git branch -a
```

### **Create Worktree for Branch**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
git worktree add ../worktrees/branch-name branch-name
```

### **Remove Worktree**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
git worktree remove ../worktrees/branch-name
```

### **Fix Detached HEAD**
```bash
git checkout -b new-branch-name
```

---

## 🚨 **COMMON ISSUES & SOLUTIONS**

### **Issue 1: "Worktree is locked"**

**Error:** `fatal: 'path' is already a working tree`

**Solution:**
```bash
# Remove lock file
rm /Users/fahadkiani/Desktop/development/crispr-assistant-main/.git/worktrees/nkw/locked
```

### **Issue 2: "Cannot remove worktree with uncommitted changes"**

**Error:** `fatal: 'path' contains modified or untracked files`

**Solution:**
```bash
# Force remove (discards changes)
git worktree remove --force ../worktrees/branch-name

# Or commit/stash first
git add . && git commit -m "Save work"
git worktree remove ../worktrees/branch-name
```

### **Issue 3: "Branch already checked out in another worktree"**

**Error:** `fatal: 'branch-name' is already checked out`

**Solution:**
- Switch to different branch in that worktree first
- Or use `git worktree add --force` (not recommended)

---

## ✅ **SUMMARY**

**Your Setup:**
- ✅ **Main repo:** `/Users/fahadkiani/Desktop/development/crispr-assistant-main` (on `main`)
- ✅ **Worktree 1:** `/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nkw` (detached HEAD)
- ✅ **Worktree 2:** `/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/rxk` (detached HEAD)

**Key Concepts:**
- Worktrees = Multiple working directories for same repo
- Each worktree can be on different branch/commit
- All worktrees share same Git history and remotes
- Cursor can open each worktree as separate workspace

**Next Steps:**
1. Fix detached HEAD in `nkw` worktree (create branch or switch to existing)
2. Commit your current changes
3. Decide: Keep worktrees or merge back to main?

---

**Questions?** Check Git documentation: `git help worktree`

