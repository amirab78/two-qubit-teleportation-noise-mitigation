# Publishing on GitHub and archiving on Zenodo (DOI)

The repository is ready as a local git repository with a clean history and the tag `v1.0.0`; nothing has been pushed.
These steps publish it on GitHub and connect it to Zenodo, which archives every GitHub release and gives it a DOI.
The whole archive, code and data together, is about 4 MB, well within the limits of both services.

## 1. Before pushing

The author, ORCID, affiliation, GitHub user and release date are filled in. The only placeholder left is `[DOI]`,
which is filled in step 6. If you change any metadata before the first push, commit it and move the tag (it has
never been pushed, so this is safe):

```bash
git commit -am "Update metadata"
```

```bash
git tag -d v1.0.0
```

```bash
git tag -a v1.0.0 -m "Version 1.0.0"
```

## 2. Create the GitHub repository

On github.com: **New repository** → name `teleportation-counts-correction`, visibility **Public**. Do **not** add a
README, licence or `.gitignore` (the repository already has them).

## 3. Push

```bash
git remote add origin https://github.com/<your-user>/teleportation-counts-correction.git
```

```bash
git push -u origin main
```

Then open the **Actions** tab and wait for the CI run to pass (tests, the hardware reproduction and all notebooks, on
Python 3.11 and 3.13; about 10 minutes).

## 4. Connect Zenodo

1. Log in to [zenodo.org](https://zenodo.org) with your GitHub account (Log in → GitHub).
2. Open your name (top right) → **GitHub**, find `teleportation-counts-correction` in the list (use **Sync now** if
   it is missing) and switch it **On**.

Zenodo takes the record's metadata (title, authors, description, licence, keywords) from `.zenodo.json`.

## 5. Make the release

Push the tag, then create the release from it:

```bash
git push origin v1.0.0
```

On GitHub: **Releases** → **Draft a new release** → choose the tag `v1.0.0` → title `v1.0.0` → paste the text of
`CHANGELOG.md` → **Publish release**. Within a few minutes Zenodo archives the release and shows it in your Zenodo
**GitHub** page with two DOIs:

* the **version DOI** (this exact release, v1.0.0): cite this one, because the results
  were produced with exactly this code and data;
* the **concept DOI** (all versions; it always resolves to the latest one).

## 6. Put the DOI into the repository

Replace `[DOI]` with the version DOI (`10.5281/zenodo.<number>`) in `README.md` (badge and citation) and
`CITATION.cff`. Add the date and the DOI to `CHANGELOG.md`, then commit and push:

```bash
git commit -am "Add the Zenodo DOI"
```

```bash
git push
```

The archived v1.0.0 itself keeps the placeholder (a release cannot contain its own DOI). To have an archived copy that
shows the DOI, make a release `v1.0.1` the same way (step 5). Zenodo then gives it a new version DOI under the same
concept DOI. In publications, keep citing v1.0.0.

## 7. Cite it

Use the BibTeX entry of the README with the version DOI, or the **Cite this repository** button that GitHub shows from
`CITATION.cff`. Zenodo's record page also offers citation formats.

## Notes

* Every commit records the name and e-mail address of your git configuration, and GitHub shows them publicly. To
  publish with GitHub's private no-reply address instead (GitHub → Settings → Emails), set it before pushing and
  rewrite the author of the local commits:
  `git config user.email "<id>+<user>@users.noreply.github.com"`, then
  `git rebase --root --exec "git commit --amend --no-edit --reset-author"`, and recreate the tag as in step 1.
* The hardware data (`data/`) are part of the same Zenodo record as the code. A separate data deposit is not needed.
* The screenshots of the IBM Quantum platform kept with the research project are not part of this repository; leave
  them out of any upload.
* The tests check that no file contains credentials, e-mail addresses, IBM instance identifiers or local paths. Run
  `pytest tests/test_repo.py` again after any change to the metadata.
* Alternative without the GitHub integration: on Zenodo, **New upload** → reserve a DOI → put it into the files →
  upload a zip of the repository at the tag.
