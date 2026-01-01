#!/bin/bash
# Delete repositories from GITHUB_REPOSITORIES.md lines 43-72 and all forks

USERNAME="fjkiani"
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Repos from lines 43-72
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

echo -e "${GREEN}✅ Found $(echo "$FORKS" | wc -l | tr -d ' ') forks${NC}"
echo -e "${GREEN}✅ Found ${#REPOS_TO_DELETE[@]} repos from lines 43-72${NC}"
echo ""
echo -e "${RED}⚠️  WARNING: This will delete $(($(echo "$FORKS" | wc -l | tr -d ' ') + ${#REPOS_TO_DELETE[@]})) repositories!${NC}"
read -p "Type 'DELETE Afirm: " CONFIRMATION

if [ "$CONFIRMATION" != "DELETE ALL" ]; then
    echo -e "${YELLOW}❌ Deletion cancelled.${NC}"
    exit 1
fi

# Delete repos from lines 43-72
echo ""
echo -e "${YELLOW}🗑️  Deleting repos from lines 43-72...${NC}"
for repo in "${REPOS_TO_DELETE[@]}"; do
    echo -e "${YELLOW}Deleting: $repo${NC}"
    gh repo delete "$USERNAME/$repo" --yes 2>&1 | grep -v "Warning:" || echo -e "${RED}Failed or already deleted: $repo${NC}"
done

# Delete all forks
echo ""
echo -e "${YELLOW}🗑️  Deleting all forks...${NC}"
echo "$FORKS" | while read -r fork; do
    if [ -n "$fork" ]; then
        echo -e "${YELLOW}Deleting fork: $fork${NC}"
        gh repo delete "$USERNAME/$fork" --yes 2>&1 | grep -v "Warning:" || echo -e "${RED}Failed or already deleted: $fork${NC}"
    fi
done

echo ""
echo -e "${GREEN}✅ Deletion process complete!${NC}"
