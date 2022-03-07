import argparse
import requests
from requests.auth import HTTPBasicAuth

# use authentication, if you want to get a lot of results. Github limits output for unauthenticated users. See https://docs.github.com/en/rest/overview/other-authentication-methods#basic-authentication
parser = argparse.ArgumentParser()
parser.add_argument('-u', '--user', type=str)
parser.add_argument('-t', '--token', type=str)
args = parser.parse_args()

organization = "precice"

commit_authors = set()
skipped_commits = []
number_commits = 0

# Collect all repos under organization
repos = set()
response = requests.get("https://api.github.com/users/{}/repos?per_page=100".format(organization))

for repo in response.json():
    repos.add(repo["name"])

print("Parsing {} repos: {}".format(len(repos), repos))

# Traverse all repos
for repo in repos:
    print("Parsing {}...".format(repo))
    pageno = 1
    api_url = "https://api.github.com/repos/{}/{}/commits?per_page=100&page={}".format(organization, repo, pageno)
    
    # Use authentication, if provided
    if args.user:
        response = requests.get(api_url, auth=HTTPBasicAuth(args.user,args.token))
    else:
        response = requests.get(api_url)
    
    # While page has entries, continue
    while len(response.json()) > 0:
        print("On page {} got {} hits.".format(pageno, len(response.json())))
        pageno += 1
        api_url = "https://api.github.com/repos/{}/{}/commits?per_page=100&page={}".format(organization, repo, pageno)
        response = requests.get(api_url, auth=HTTPBasicAuth(args.user,args.token))
        # Traverse all commits on page
        for commit in response.json():
            number_commits += 1
            try:
                # extract name of author
                author = commit['author']['login']
                commit_authors.add(author)
            except TypeError:
                # need to skip many commits. Not sure why.
                skipped_commits.append(commit)

print("Parsed {} commits in total.".format(number_commits))
print("Found {} unique authors: {}".format(len(commit_authors), commit_authors))
print("Skipped {} commits due to parsing errors".format(len(skipped_commits)))
