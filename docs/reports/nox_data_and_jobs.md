# Nox Data and Jobs Inventory

This document summarizes all database models, migrations, SQL/PLpgSQL files, and background job/worker queue systems in the `nox` repository.  
Sources: migration scripts, milestone reports, technical docs, and extracted SQL.

---

## Database Models

### Core Tables
- **users**
  - Fields: id (UUID, PK), email, role, is_active, quota_files, quota_cpu_seconds, quota_memory_mb, quota_req_hour, quota_req_day, quota_storage_mb, quota_files_max
  - Relationships: Linked to user_usage, quota_violations, batch_job_collections
- **user_usage**
  - Tracks real-time resource usage per user
  - Foreign key: user_id → users(id)
- **quota_violations**
  - Records quota violations (JSONB details)
  - Foreign key: user_id → users(id)
- **audit_sessions**
  - Session management (UUID, user_id, client_ip, user_agent, timestamps)
  - 9 indexes for performance
- **audit_actions**
  - Detailed action logs (session_id, action, resource, timestamps)
  - 11 indexes
- **audit_admin_actions**
  - Admin action tracking (action, resource_type, resource_id)
  - 3 indexes
- **audit_daily_summaries**
  - Daily rollups, automated by triggers
  - 2 indexes + trigger

### OAuth2 and Security
- **oauth2_tokens**
  - Token storage, refresh logic
  - Indexed for fast lookup
- **oauth2_profiles**
  - Comprehensive user profiles
- **oauth2_login_sessions**
  - Tracks login events, audit integration
- **oauth2_token_refreshes**
  - Refresh and revocation events
- **oauth2_system_config**
  - Configuration settings

### Jobs, Batching, Metrics
- **computational_jobs**
  - Fields: id (UUID, PK), job_version, error_recovery_attempts, parent_job_id (FK → computational_jobs.id)
  - Relationships: batch_job_collections
- **batch_job_collections**
  - Fields: id, user_id (FK → users.id), collection_name, total_jobs, completed_jobs, failed_jobs, status, created_at, started_at, completed_at
- **webvitals_metrics**
  - Web performance metrics, for operational monitoring
- **ai_security_audit**
  - AI security audit logs
- **websocket_connections**
  - Tracks live WebSocket sessions
- **system_configuration**
  - Stores system-wide config parameters
- **admin_audit_trail**
  - Tracks admin actions and resources
- **rate_limit_buckets**
  - For API rate limiting

---

## Migrations

### Migration File Order
1. **v7.x → v8.0.0**  
   - File: `migrations/v8_to_v8.0.0.sql`  
   - Rollback and validation scripts present

2. **User Quotas & Usage (M5.1)**
   - Adds quota columns to users, user_usage and quota_violations tables

3. **M6 Audit Logging**
   - Adds audit_sessions, audit_actions, audit_admin_actions, audit_daily_summaries

4. **M7 OAuth2 Integration**
   - Adds oauth2_tokens, oauth2_profiles, oauth2_login_sessions, oauth2_token_refreshes, oauth2_system_config

5. **WebVitals & AI Security (v8.0.0)**
   - Adds webvitals_metrics, ai_security_audit, websocket_connections

---

## Stored Procedures, Functions, Triggers

- **audit_daily_summaries**: Daily rollup via PostgreSQL triggers
- **oauth2_system_config**: 3 utility functions for config management
- **admin_audit_trail**: 3 automated triggers for admin action logging

---

## Background Jobs & Schedulers

### Job Processing
- **computational_jobs**: Core job tracking table
- **batch_job_collections**: Groups related jobs for batch processing
- **Background processing**: Likely uses Python asyncio or Celery (not explicitly found in search results)

### Queue Systems
- **Rate limiting**: rate_limit_buckets table for API throttling
- **WebSocket management**: websocket_connections for real-time updates
- **Job recovery**: error_recovery_attempts field in computational_jobs

### Monitoring & Metrics
- **webvitals_metrics**: Frontend performance monitoring
- **ai_security_audit**: Security event tracking
- **audit_daily_summaries**: Automated daily rollups via triggers

---

## Database Configuration

- **Primary Database**: PostgreSQL (inferred from PLpgSQL usage)
- **Indexing Strategy**: Heavy indexing on audit tables (9-11 indexes per table)
- **Automation**: Triggers for daily summaries and audit trail
- **Security**: OAuth2 integration with token management

---

**Note**: Background job processing details may be in files not captured in search results. Check for Celery, RQ, or custom job runner configurations in the codebase.
