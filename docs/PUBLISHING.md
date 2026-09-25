# Manual publication only

No Git remote operation was performed. `GITHUB_UPLOAD_READY` names the requested code-only candidate directory, not a full-reproduction certification. Read PORTABILITY.md before publication and describe limitations in the release.

After the owner authorizes publishing, clone the existing repository separately and create a new branch. Back up its old tree; preserve `.git`, valid existing licences and all v1 history. Copy only the contents of this candidate directory. Do not copy LOCAL_REVIEW_ONLY, the execution task, input ZIPs or the delivery ZIP.

If the old root `code.zip` is tracked, only after its contents are migrated/backed up remove that single file from the new index using `git rm --cached code.zip`. Do not clear the entire index or rewrite history. Review the staged diff for data/results/credentials. Push and create a new release only with explicit authorization. Update manuscript code availability only after an actual remote version is reachable; this task did not change the manuscript.
