# `CLAUDE.md` - nf-metabarcoding


## The Golden Rule

When unsure about implementation details, ALWAYS ask the developer.


## Project Context

nf-metabarcoding is an open-source nextflow workflow for microbiome
analysis.


### Guidelines:

- **always ask a human for instruction** if you think a comment should
  be modified or removed, or if you think a non-target test file
  should be modified
- in rare occasions the specified behaviour (see SPECIFICATIONS.md and
  the actual behaviour of the workflow may differ. If you identify
  such a situation, it is very important that you describe the case
  and **ask** for a human review
- if the README or SPECIFICATIONS lack precision regarding the actual
  behaviour, you can suggest improvements to these files
- if this file (CLAUDE.md) is unclear, or if you think a change would
  yield better refactoring results, feel free to suggest improvements
- test-driven development (write tests first, each specification must
  be covered by tests)


### Refactoring instructions

In this environment, you have access to `git`, and you have
branch-creation rights:
- create a new git branch, stemming from the `dev` branch, and using
  the naming pattern `tmp_$(date +%Y%m%d%H%M%S)`
- checkout that new branch
- commit often, it is ok to have small changes in each commit, it
  makes bug-hunting easier
- add "Co-Authored-By" trailers to commit messages, for my co-worker
  Florian Filloux (github handle: ffillouxdev)
- do not add yourself (Claude) as co-author in commit messages (sorry)
- when done, ask for a human review. Humans are in charge of merging
  your work into the dev branch
- check modified bash files with `shellcheck` and fix reported
  issues.
- check modified python files with `flake` and fix reported issues.


### Bash conventions

- declare one variable per line — do not group declarations such as
  `local a="x" b="y"` or `local x; x="..."` on a single line
- mark variables read-only when possible:
  - top-level constants use `readonly NAME=value` (or `readonly -a`
    for arrays)
  - single-shot function-scope variables use `local -r`
  - inside loop bodies, use plain `local` — `local -r` fails on the
    second iteration because the variable is already read-only
- when assigning the result of a command substitution to a local
  variable, split it onto two lines so the command's exit status is
  not masked by `local`'s (see shellcheck SC2155):

  ```bash
  local q
  q="$(some_command)"
  ```


## What AI Must always Do

AI must always leave these untouched:
1. **comments** - They're there for a reason
2. **already existing tests** - They've been reviewed by a human
3. **CLAUDE.md**

Remember: We prefer maintainability, correctness, and readability over
cleverness. When in doubt, choose the boring solution or ask for help.
