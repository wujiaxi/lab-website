# 🔬 Lab Website - Local Development

A Docker-based local WordPress development environment for your lab website.

## 🚀 Quick Start

### 1. Start the containers
```bash
docker-compose up -d
```

### 2. Access your site
- **WordPress**: http://localhost:8080
- **phpMyAdmin**: http://localhost:8081 (for database management)

### 3. Complete WordPress Setup
1. Open http://localhost:8080 in your browser
2. Select your language
3. Fill in site information:
   - Site Title: Your Lab Name
   - Username: admin
   - Password: (choose a strong password)
   - Email: your@email.com
4. Click "Install WordPress"

## 📁 Project Structure

```
lab-website/
├── docker-compose.yml      # Docker configuration
├── README.md               # This file
└── wp-content/
    ├── themes/             # Custom themes go here
    ├── plugins/            # Custom plugins go here
    └── uploads/            # Media uploads
```

## 🛠️ Useful Commands

```bash
# Start containers
docker-compose up -d

# Stop containers
docker-compose down

# View logs
docker-compose logs -f

# Restart WordPress
docker-compose restart wordpress

# Access WordPress container shell
docker exec -it lab-website-wp bash

# Backup database
docker exec lab-website-db mysqldump -u wordpress -pwordpress_password wordpress > backup.sql
```

## 🎨 Recommended Next Steps

1. **Install a Theme**: Go to Appearance → Themes → Add New
2. **Install Plugins**: Go to Plugins → Add New
   - Elementor (page builder)
   - Jenga (publications management)
   - WPDataTables (dynamic data display)
3. **Create Pages**: Pages → Add New

## 🔧 Troubleshooting

### Port already in use?
Change the port in `docker-compose.yml`:
```yaml
ports:
  - "8888:80"  # Change 8080 to another port
```

### Reset everything?
```bash
docker-compose down -v  # This removes all data!
docker-compose up -d
```

## 📊 Database Credentials

| Setting | Value |
|---------|-------|
| Host | db |
| Database | wordpress |
| User | wordpress |
| Password | wordpress_password |
