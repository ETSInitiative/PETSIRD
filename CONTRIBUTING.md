Contributing
============

(This document largely follows the [SIRF Contribution guidelines](https://github.com/SyneRBI/SIRF/blob/master/CONTRIBUTING.md).)

Please help us by finding problems, discussing on the mailing lists, contributing documentation,
bug fixes or even features. Below are some brief guidelines.

## Reporting a problem

Please use our [issue-tracker].

## Submitting a patch

For contributing any code or documentation that is non-trivial, we require a
signed Contributor License Agreement, stating clearly that your
conributions are licensed appropriately. This will normally need to be signed by your
employer/university, unless you own your own copyright.
You will have to do this only once. Please contact us for more information.

Please keep a patch focused on a single issue/feature. This is important to keep our history clean,
but will also help reviewing things and therefore speed-up acceptance.

### Process

This is our recommended process. If it sounds too daunting, ask for help.

1. Create a new issue (see above). State that you will contribute a fix if you intend to do so.
2. Create a [fork](https://help.github.com/articles/fork-a-repo) on github and work from there.
3. Create a branch in your fork with a descriptive name and put your fixes there. If your fix is
simple you could do it on github by editing a file, otherwise clone your project (or add a remote
to your current git clone) and work as usual.
4. If your change is important, add yourself to `NOTICE.txt`, describe it in `CHANGES.md` (if it exists)
and other documentation files.
5. Use [well-formed commit messages](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
for each change (in particular with a single "subject" line
followed by an empty line and then more details).
6. Push the commits to your fork and submit a [pull request (PR)](https://help.github.com/articles/creating-a-pull-request)
(enable changes by project admins.) Be prepared to add further commits to your branch after discussion.
In the description of the PR, add a statement about which Issue this applies to
using [a phrase such that github auto-closes the issue when merged to master](https://help.github.com/articles/closing-issues-using-keywords/).
7. Be prepared to add further commits to your branch after discussion.
Please by mindful about the resources used by our Continuous Integration (CI) workflows:
  - Group your commits and only push once your code compiles and tests succeed on your machine
  - Use specific keywords in the body of the last commit that you push to prevent CI being run:
     - `[ci skip]` skips all CI runs (e.g. when you only change documentation, or when your update isn't ready yet)
     - `[actions skip]` does not run GitHub Actions, see [here](https://github.blog/changelog/2021-02-08-github-actions-skip-pull-request-and-push-workflows-with-skip-ci/).
8. After acceptance of your PR, go home with a nice warm feeling.

Suggested reading:
https://help.github.com/articles/fork-a-repo/, https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project or https://guides.github.com/activities/forking/.

### A note on copyright dates and notices (and licenses)

We recommend adding a brief copyright and license header to new files, e.g. in C++
```c++
/*
  Copyright (C) 2023 University College London

  SPDX-License-Identifier: Apache-2.0
*/
```
If you modify an existing file, you need to make sure the copyright header is up-to-date for your changes
(unless it's a trivial change).

If you copied code from somewhere, you need to preserve its copyright date/notice. If you copied external code,
you need to make sure its license is compatible with the Apache 2.0 license, and indicate clearly what the license
of the copied code is (and follow its terms of course).

## Project rules

- Only one official, stable, up-to-date branch: **main**
    + Essentially "latest stable beta version with no known bugs
      since the last official release version"
    + Never knowingly add a bug to **main**
- Any work-in-progress commits should be in their own branches.
- GitHub assigns a unique number to each issue, c.f. the [issue-tracker].
- A pull request (PR) is an issue with an associated branch,
  c.f. [pull-requests]. Even for "internal" development, we prefer a PR for
  a branch to allow review and discussion.
- Branches and PRs are kept small (ideally one 'feature' only) and branch from **main**,
  not from another branch, unless required. This allows
  commenting/improving/merging this branch/PR
  independent of other developments.
- Contributions of new features should also update documentation and release notes. After version 1.0,
  this documentation needs to state something like "introduced after version 1.xxx".
- We prefer issues to be opened via [github][issue-tracker] due to the following reasons:
    + Ensures issues will never get lost in emails
        * Facilitates issue status tracking
    + Allows focused comments/discussion
        * Easy cross-referencing of related issues, PRs, and commits

```

[issue-tracker]: https://github.com/ETSInitiative/PRDdefinition/issues
[submodules]: https://git-scm.com/docs/gitsubmodules
[pull-requests]: https://github.com/ETSInitiative/PRDdefinition/pulls
