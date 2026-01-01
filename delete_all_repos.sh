#!/bin/bash
# Delete repositories from GITHUB_REPOSITORIES.md lines 43-72 and all forks

USERNAME="fjkiani"
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Repos from lines 43-72 (excluding crispro - we're using that one)
REPOS_TO_DELETE=(
    "crm-deployment"
    "deepfake-bitmind"
    "defiKSA"
    "awesome-cursorrules"
    "claude-code"
    "life-sciences"
    "genomics-mcp"
    "data-room"
    "oncology-backend-v2"
    "jedi-v2"
    "fahad-kiani"
    "sunvic-2"
    "SunVic"
    "lotto-machine"
)

# Get all forks
echo -e "${YELLOW}📋 Fetching all forks...${NC}"
FORKS=$(gh repo list $USERNAME --limit 500 --json name,isFork --jq '.[] | select(.isFork == true) | .name')

FORK_COUNT=$(echo "$FORKS" | grep -v '^$' | wc -l | tr -d ' ')
REPO_COUNT=${#REPOS_TO_DELETE[@]}
TOTAL=$((FORK_COUNT + REPO_COUNT))

echo -e "${GREEN}✅ Found $FORK_COUNT forks${NC}"
echo -e "${GREEN}✅ Found $REPO_COUNT repos from lines 43-72${NC}"
echo -e "${RED}⚠️  Wis will delete $TOTAL repositories${NC}"
echo ""
read -p "Type 'DELETE ALL' to confirm: " CONFIRMATION

if [ "$CONFIRMATION" != "DELETE ALL" ]; then
    echo -e "${YELLOW}Deletion cancelled.${NC}"
    exit 0
fi

# Delete repos from file
echo ""
echo -e "${YELLOW}🗑️  Deleting repos from lines 43-72...${NC}"
for repo in "${REPOS_TO_DELETE[@]}"; do
    echo -e "${YELLOW}Deleting: $USERNAME/$repo${NC}"
    gh repo delete "$USERNAME/$repo" --yes 2>&1 | grep -v "Warning:" || echo -e "${RED}Failed or already deleted: $repo${NC}"
done

# Delete all forks
echo ""
echo -e "${YELLOW}🗑️  Deleting all forks...${NC}"
DELETED=0
FAILED=0
while IFS= read -r fork; do
    if [ -n "$fork" ]; then
        echo -e "${YELLOW}Deleting fork: $USERNAME/$fork${NC}"
        if gh repo delete "$USERNAME/$fork" --yes 2>&1 | grep -v "Warning:"; then
            DELETED=$((DELETED + 1))
        else
            FAILED=$((FAILED + 1))
            echo -e "${RED}Failed or already deleted: $fork${NC}"
        fi
    fi
done <<< "$Fcho ""
echo -e "${GREEN}✅ Deletion complete!${NC}"
echo -e "${GREEN}   Deleted: $DELETED forks${NC}"
echo -e "${RED}   Failed: $FAILED forks${NC}"
