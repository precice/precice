#!/bin/bash

# Checks pull request to match certain style and post in the
# pull request thread if this style is not matched

# Checks code style with clang-tidy
check_code_style()
{
  not_tidy='   '
  # check results with clang-tidy
  cd $TRAVIS_BUILD_DIR
  cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON $PRECICE_ROOT > /dev/null
  for file in ${files_changed}; do
    # clang-tidy onyl reports on files that are compiled
    if [ ${file##*.} == 'cpp' ]; then
      clang-tidy -p $TRAVIS_BUILD_DIR $file
      if [ "$?" -eq 0 ]; then
        not_tidy+="\n    * \`$file\` "
        pr_invalid=1
      fi
    fi
  done

  if [ -n "$not_tidy" ]; then
    bot_msg+="\n* Clang-tidy complained on the following files: "
    bot_msg+="${not_tidy}"
    bot_msg+="\n Consider running clang-tidy on them ( or browsing [Travis log](https://api.travis-ci.org/v3/job/${TRAVIS_JOB_ID}/log.txt ) ) and fixing reported issues"
  fi
}


# Checks code formatting with clang-format
check_code_format()
{
  not_formatted='   '
  # Get diff from this commit and filter filetypes that are needed
  if [ -n "${files_changed}" ]; then
    for file in ${files_changed}; do
      clang-format -style=file -output-replacements-xml $file  | grep -c "<replacement " > /dev/null
      if [ "$?" -eq 0 ]; then
        not_formatted+="\n    * \`$file\` "
        pr_invalid=1
      fi
    done
  fi

  if [ -n "$not_formatted" ]; then
    bot_msg+="\n* Your code formatting did not follow our clang-format style in following files: "
    bot_msg+="$not_formatted"
  fi
}

# Checks if something was added to the changelog, if certain threshold
# of lines  was changed
check_changelog()
{
  # perform merge and get the number of lines changed
  lines_changed=$(git log | sed -n '2p' | awk '{print $2, $3}' | xargs git diff --numstat | awk '{ sum+= $1 + $2 ; } END { print sum; }')
  if [ $lines_changed -gt 100 ]; then
    git log | sed -n '2p' | awk '{print $2, $3}' | xargs git diff --numstat | grep 'CHANGELOG.md'
    if [ "$?" -eq 1 ]; then
      pr_invalid=1
      bot_msg+="\n* It seems like you  forgot to update \`CHANGELOG.md\`"
    fi
  fi
}

# Get list of files that were changed
files_changed=$( git log | sed -n '2p' | awk '{print $2, $3}' | xargs git diff --name-only | grep '.cpp\|.hpp' | xargs )
pr_invalid=0
# generic bot message start
bot_msg="Thank you for your contribution.\n"

# perform necessary checks
check_changelog
check_code_format
check_code_style

# Send message to github API to post in the PR thread if we failed
if [[ "$pr_invalid" -eq 1 ]]; then
   curl -s -H "Authorization: token $TRAVIS_ACCESS_TOKEN" -X POST -d "{\"body\": \"$bot_msg\"}" \
     "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/issues/${TRAVIS_PULL_REQUEST}/comments" > /dev/null
   exit 1
fi
