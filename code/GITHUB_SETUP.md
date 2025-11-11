# GitHub Setup Guide

This guide will walk you through setting up your Wastewater analysis repository on GitHub.

## Prerequisites

- GitHub account (free at https://github.com/signup)
- Git installed locally (https://git-scm.com/downloads)
- The Wastewater analysis code and data ready locally

## Step 1: Create GitHub Repository

1. Go to https://github.com/new
2. Fill in repository details:
   - **Repository name**: `Wastewater`
   - **Description**: "FoodSeq-FLOW Dietary DNA Analysis from Wastewater"
   - **Visibility**: Public (recommended for research transparency)
   - **Initialize**: Leave unchecked (we'll push existing code)
3. Click **Create repository**
4. Copy the repository URL (should be: `https://github.com/YOUR_USERNAME/Wastewater.git`)

## Step 2: Configure Git Locally

Open terminal/command prompt and run:

```bash
# Set your GitHub credentials (one time)
git config --global user.name "Your Name"
git config --global user.email "your.email@github.com"

# If you want to avoid password prompts, set up SSH:
# https://docs.github.com/en/authentication/connecting-to-github-with-ssh
```

## Step 3: Prepare Local Repository

Navigate to your project directory:

```bash
cd "C:/Users/rache/Box/project_davidlab/LAD_LAB_Personnel/Rachel_Q/Code and Data"

# Initialize git (if not already done)
git init

# Check status
git status
```

## Step 4: Add Files to Git

```bash
# Add all files (respects .gitignore)
git add .

# Or add specific directories
git add NCWW_ms_code/
git add "Data for code optimization_Do not submit/"
git add ANALYSIS_REPRODUCTION_GUIDE.md

# Check what will be committed
git status
```

## Step 5: Create Initial Commit

```bash
git commit -m "Initial commit: Complete FoodSeq-FLOW analysis pipeline

- Added run_full_analysis.R script for reproducible analysis
- Added comprehensive documentation (README, QUICK_START, ANALYSIS_REPRODUCTION_GUIDE)
- Added all data files (phyloseq objects, metadata, taxonomy)
- Added original RMarkdown analysis code
- Reproducible analysis of wastewater dietary DNA from 183 NC samples"
```

## Step 6: Connect to GitHub

```bash
# Add GitHub repository as remote (replace YOUR_USERNAME)
git remote add origin https://github.com/YOUR_USERNAME/Wastewater.git

# Verify connection
git remote -v
# Should show:
# origin  https://github.com/YOUR_USERNAME/Wastewater.git (fetch)
# origin  https://github.com/YOUR_USERNAME/Wastewater.git (push)

# Rename branch to main (GitHub default)
git branch -M main

# Push to GitHub
git push -u origin main
```

## Step 7: Verify on GitHub

1. Go to https://github.com/YOUR_USERNAME/Wastewater
2. You should see:
   - All code files
   - README.md displayed as project description
   - File tree showing structure
   - Git history showing your commit

## Step 8: Add GitHub Badges (Optional)

Edit README.md and add badges after the title:

```markdown
# FoodSeq-FLOW Wastewater Dietary DNA Analysis

[![GitHub stars](https://img.shields.io/github/stars/YOUR_USERNAME/Wastewater?style=flat-square)](https://github.com/YOUR_USERNAME/Wastewater)
[![GitHub license](https://img.shields.io/github/license/YOUR_USERNAME/Wastewater?style=flat-square)](LICENSE)
```

Then commit:
```bash
git add README.md
git commit -m "Add GitHub badges to README"
git push
```

## Pushing Updated Results

After running the analysis locally and generating figures:

```bash
# Add new figures and results
git add NCWW_ms_code/figures/
git add NCWW_ms_code/results/

# Commit with description
git commit -m "Add analysis results: Generated figures and statistics tables

- 8 publication-quality figures (PNG, 300 DPI)
- Summary statistics and diversity metrics tables
- All figures show expected results from analysis"

# Push to GitHub
git push
```

## Troubleshooting

### Issue: "fatal: not a git repository"
**Solution**: Run `git init` in your project directory first

### Issue: "Permission denied (publickey)"
**Solution**: Set up SSH key:
1. Generate key: `ssh-keygen -t ed25519 -C "your_email@example.com"`
2. Add to GitHub: Settings → SSH and GPG keys → New SSH key
3. Test: `ssh -T git@github.com`

### Issue: "authentication failed"
**Solution A**: Use token authentication:
```bash
git remote set-url origin https://YOUR_TOKEN@github.com/YOUR_USERNAME/Wastewater.git
```
Get token at: Settings → Developer settings → Personal access tokens

**Solution B**: Use GitHub Desktop (GUI):
- Download: https://desktop.github.com/
- Much easier for beginners

### Issue: Large files rejected by GitHub
**Solution**: Use Git LFS for large files:
```bash
git lfs install
git lfs track "*.rds"
git add .gitattributes
git commit -m "Add LFS tracking for large .rds files"
git push
```

## Repository Structure on GitHub

After pushing, your repository will look like:

```
Wastewater/
├── README.md                          ← Main project description
├── QUICK_START.md                     ← 5-minute setup guide
├── GITHUB_SETUP.md                    ← This file
├── ANALYSIS_REPRODUCTION_GUIDE.md     ← Detailed methods
├── .gitignore                         ← Files to exclude
├── NCWW_ms_code/
│   ├── run_full_analysis.R            ← Main analysis script
│   ├── Rcode.Rmd                      ← Original code
│   ├── Rcode_reproduced.R             ← Alternative version
│   ├── data/                          ← All data files
│   ├── figures/                       ← Generated plots
│   └── results/                       ← Generated tables
└── Data for code optimization_Do not submit/
    └── [Raw data files]
```

## Collaborate with Others

To allow others to contribute:

1. Go to repository Settings → Collaborators → Add people
2. Or create a pull request workflow:
   - Others fork the repository
   - Make changes in their fork
   - Submit pull request for review
   - You merge approved changes

## Additional GitHub Features

### Enable GitHub Pages (Host Results Online)
1. Settings → Pages → Source → main branch
2. Your README becomes a website at: `https://YOUR_USERNAME.github.io/Wastewater/`

### Create Releases
```bash
git tag -a v1.0 -m "Initial analysis release"
git push origin v1.0
```
Then create release on GitHub with:
- Description of changes
- Link to manuscript
- Download links for figures

### Add Discussion Forum
1. Settings → Features → Check "Discussions"
2. Others can ask questions about the analysis

### Enable GitHub Actions (CI/CD)
Create `.github/workflows/test.yml`:
```yaml
name: Test Analysis

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: r-lib/actions/setup-r@v2
      - run: Rscript -e "source('run_full_analysis.R')"
```

## Citing Your Repository

Add to your README:

```bibtex
@software{Wastewater2025,
  author = {Your Name},
  title = {Wastewater Dietary DNA Analysis},
  url = {https://github.com/YOUR_USERNAME/Wastewater},
  year = {2025}
}
```

## Common Git Commands

```bash
# View commit history
git log --oneline

# View changes not yet staged
git diff

# View staged changes
git diff --staged

# Undo last commit (keep changes)
git reset HEAD~1

# Add changes interactively
git add -p

# View who changed each line
git blame filename.R

# Create a branch for new feature
git checkout -b feature/new-analysis

# Switch branches
git checkout main

# Delete a branch
git branch -d feature/new-analysis
```

## Next Steps

1. ✓ Create GitHub repository
2. ✓ Push code and documentation
3. **Run analysis locally** to generate figures
4. Push figures and results to GitHub
5. Share repository URL with collaborators
6. Enable GitHub Pages for public access (optional)
7. Monitor pull requests and issues

## Resources

- GitHub Docs: https://docs.github.com
- Git Cheat Sheet: https://github.github.com/training-kit/downloads/github-git-cheat-sheet.pdf
- Markdown Guide: https://www.markdownguide.org/
- R + Git Integration: https://happygitwithr.com/

---

**Need help?**
- GitHub Support: https://support.github.com
- Contact: YOUR_EMAIL@example.com
