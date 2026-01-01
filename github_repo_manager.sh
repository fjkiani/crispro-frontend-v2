#!/bin/bash
# GitHub Repo Manager - CLI tool for organizing fjkiani repos

USERNAME="fjkiani"
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

list_repos() {
    echo -e "${GREEN}📋 Your Repositories:${NC}"
    gh repo list $USERNAME --limit 500 --json name,isPrivate,updatedAt,isArchived,isFork \
        --jq 'sort_by(.updatedAt) | reverse | .[] | 
        "\(.name) | Private:\(.isPrivate) | Archived:\(.isArchived) | Fork:\(.isFork) | Updated:\(.updatedAt[:10])"'
}

list_forks() {
    echo -e "${YELLOW}🍴 Your Forks:${NC}"
    gh repo list $USERNAME --limit 500 --json name,isFork --jq '.[] | select(.isFork == true) | .name'
}

list_own_repos() {
    echo -e "${GREEN}⭐ Your Own Repos (not forks):${NC}"
    gh repo list $USERNAME --limit 500 --json name,isFork --jq '.[] | select(.isFork == false) | .name'
}

delete_repo() {
    if [ -z "$1" ]; then
        echo -e "${RED}❌ Error: Please provide a repo name${NC}"
        exit 1
    fi
    echo -e "${RED}⚠️  WARNING: This will permanently delete $USERNAME/$1${NC}"
    read -p "Are you sure? (yes/no): " confirm
    if [ "$confirm" = "yes" ]; then
        gh repo delete $USERNAME/$1 --yes
        echo -e "${GREEN}✅ Deleted $USERNAME/$1${NC}"
    else
        echo "Cancelled."
    fi
}

archive_repo() {
    if [ -z "$1" ]; then
        echo -e "${RED}❌ Error: Please provide a repo name${NC}"
        exit 1
    fi
    gh repo archive $USERNAME/$1
    echo -e "${GREEN}✅ Archived $USERNAME/$1${NC}"
}

create_repo() {
    if [ -z "$1" ]; then
        echo -e "${RED}❌ Error: Please provide a repo name${NC}"
        exit 1
    fi
    if [ "$2" = "--private" ]; then
        gh repo create $1 --private
    else
        gh repo create $1 --public
    fi
    echo -e "${GREEN}✅ Created $1${NC}"
}

show_help() {
    echo "GitHub Repo Manager for $USERNAME"
    echo ""
    echo "Usage: $0 [command] [options]"
    echo ""
    echo "Commands:"
    echo "  list              List all repos"
    echo "  forks             List only forks"
    echo "  own               List only your repos (not forks)"
    echo "  delete REPO       Delete a repo (with confirmation)"
    echo "  archive REPO      Archive a repo"
    echo "  create REPO       Create a new public repo"
    echo "  create REPO --private  Create a new private repo"
    echo "  help              Show this help"
}

case "$1" in
    list) list_repos ;;
    forks) list_forks ;;
    own) list_own_repos ;;
    delete) delete_repo "$2" ;;
    archive) archive_repo "$2" ;;
    create) create_repo "$2" "$3" ;;
    help|--help|-h|"") show_help ;;
    *) echo "Unknown command: $1"; show_help; exit 1 ;;
esac

# GitHub Repo Manager - CLI tool for organizing fjkiani repos

USERNAME="fjkiani"
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

list_repos() {
    echo -e "${GREEN}📋 Your Repositories:${NC}"
    gh repo list $USERNAME --limit 500 --json name,isPrivate,updatedAt,isArchived,isFork \
        --jq 'sort_by(.updatedAt) | reverse | .[] | 
        "\(.name) | Private:\(.isPrivate) | Archived:\(.isArchived) | Fork:\(.isFork) | Updated:\(.updatedAt[:10])"'
}

list_forks() {
    echo -e "${YELLOW}🍴 Your Forks:${NC}"
    gh repo list $USERNAME --limit 500 --json name,isFork --jq '.[] | select(.isFork == true) | .name'
}

list_own_repos() {
    echo -e "${GREEN}⭐ Your Own Repos (not forks):${NC}"
    gh repo list $USERNAME --limit 500 --json name,isFork --jq '.[] | select(.isFork == false) | .name'
}

delete_repo() {
    if [ -z "$1" ]; then
        echo -e "${RED}❌ Error: Please provide a repo name${NC}"
        exit 1
    fi
    echo -e "${RED}⚠️  WARNING: This will permanently delete $USERNAME/$1${NC}"
    read -p "Are you sure? (yes/no): " confirm
    if [ "$confirm" = "yes" ]; then
        gh repo delete $USERNAME/$1 --yes
        echo -e "${GREEN}✅ Deleted $USERNAME/$1${NC}"
    else
        echo "Cancelled."
    fi
}

archive_repo() {
    if [ -z "$1" ]; then
        echo -e "${RED}❌ Error: Please provide a repo name${NC}"
        exit 1
    fi
    gh repo archive $USERNAME/$1
    echo -e "${GREEN}✅ Archived $USERNAME/$1${NC}"
}

create_repo() {
    if [ -z "$1" ]; then
        echo -e "${RED}❌ Error: Please provide a repo name${NC}"
        exit 1
    fi
    if [ "$2" = "--private" ]; then
        gh repo create $1 --private
    else
        gh repo create $1 --public
    fi
    echo -e "${GREEN}✅ Created $1${NC}"
}

show_help() {
    echo "GitHub Repo Manager for $USERNAME"
    echo ""
    echo "Usage: $0 [command] [options]"
    echo ""
    echo "Commands:"
    echo "  list              List all repos"
    echo "  forks             List only forks"
    echo "  own               List only your repos (not forks)"
    echo "  delete REPO       Delete a repo (with confirmation)"
    echo "  archive REPO      Archive a repo"
    echo "  create REPO       Create a new public repo"
    echo "  create REPO --private  Create a new private repo"
    echo "  help              Show this help"
}

case "$1" in
    list) list_repos ;;
    forks) list_forks ;;
    own) list_own_repos ;;
    delete) delete_repo "$2" ;;
    archive) archive_repo "$2" ;;
    create) create_repo "$2" "$3" ;;
    help|--help|-h|"") show_help ;;
    *) echo "Unknown command: $1"; show_help; exit 1 ;;
esac

# GitHub Repo Manager - CLI tool for organizing fjkiani repos

USERNAME="fjkiani"
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

list_repos() {
    echo -e "${GREEN}📋 Your Repositories:${NC}"
    gh repo list $USERNAME --limit 500 --json name,isPrivate,updatedAt,isArchived,isFork \
        --jq 'sort_by(.updatedAt) | reverse | .[] | 
        "\(.name) | Private:\(.isPrivate) | Archived:\(.isArchived) | Fork:\(.isFork) | Updated:\(.updatedAt[:10])"'
}

list_forks() {
    echo -e "${YELLOW}🍴 Your Forks:${NC}"
    gh repo list $USERNAME --limit 500 --json name,isFork --jq '.[] | select(.isFork == true) | .name'
}

list_own_repos() {
    echo -e "${GREEN}⭐ Your Own Repos (not forks):${NC}"
    gh repo list $USERNAME --limit 500 --json name,isFork --jq '.[] | select(.isFork == false) | .name'
}

delete_repo() {
    if [ -z "$1" ]; then
        echo -e "${RED}❌ Error: Please provide a repo name${NC}"
        exit 1
    fi
    echo -e "${RED}⚠️  WARNING: This will permanently delete $USERNAME/$1${NC}"
    read -p "Are you sure? (yes/no): " confirm
    if [ "$confirm" = "yes" ]; then
        gh repo delete $USERNAME/$1 --yes
        echo -e "${GREEN}✅ Deleted $USERNAME/$1${NC}"
    else
        echo "Cancelled."
    fi
}

archive_repo() {
    if [ -z "$1" ]; then
        echo -e "${RED}❌ Error: Please provide a repo name${NC}"
        exit 1
    fi
    gh repo archive $USERNAME/$1
    echo -e "${GREEN}✅ Archived $USERNAME/$1${NC}"
}

create_repo() {
    if [ -z "$1" ]; then
        echo -e "${RED}❌ Error: Please provide a repo name${NC}"
        exit 1
    fi
    if [ "$2" = "--private" ]; then
        gh repo create $1 --private
    else
        gh repo create $1 --public
    fi
    echo -e "${GREEN}✅ Created $1${NC}"
}

show_help() {
    echo "GitHub Repo Manager for $USERNAME"
    echo ""
    echo "Usage: $0 [command] [options]"
    echo ""
    echo "Commands:"
    echo "  list              List all repos"
    echo "  forks             List only forks"
    echo "  own               List only your repos (not forks)"
    echo "  delete REPO       Delete a repo (with confirmation)"
    echo "  archive REPO      Archive a repo"
    echo "  create REPO       Create a new public repo"
    echo "  create REPO --private  Create a new private repo"
    echo "  help              Show this help"
}

case "$1" in
    list) list_repos ;;
    forks) list_forks ;;
    own) list_own_repos ;;
    delete) delete_repo "$2" ;;
    archive) archive_repo "$2" ;;
    create) create_repo "$2" "$3" ;;
    help|--help|-h|"") show_help ;;
    *) echo "Unknown command: $1"; show_help; exit 1 ;;
esac


