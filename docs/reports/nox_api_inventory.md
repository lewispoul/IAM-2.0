# Nox API Endpoint Inventory

Below is a structured inventory of API/web endpoints, routers, middlewares, and related details as defined in the `nox` repository.

---

| HTTP Method | Path                   | Handler Function                | Source File:Line             | Request Schema       | Response Schema        | Auth/Deps/Middleware         | Notes/Docs |
|-------------|------------------------|---------------------------------|------------------------------|----------------------|------------------------|------------------------------|------------|
| GET         | `/me`                  | get_current_user_info           | auth/routes.py:73            | (Depends current_user) | UserOut               | Depends(get_current_user)    | User info  |
| GET         | `/users`               | list_users                      | auth/routes.py:78            | limit, offset        | List[UserOut]          | Depends(require_role(ADMIN)) | Admin only |
| GET         | `/users/{user_id}`     | get_user                        | auth/routes.py:105           | user_id              | UserOut                | Depends(require_role(ADMIN)) | Admin only |
| PUT         | `/users/{user_id}`     | update_user                     | auth/routes.py:117           | user_id, UserUpdate  | UserOut                | Depends(require_role(ADMIN)) | Admin only |
| DELETE      | `/users/{user_id}`     | delete_user                     | auth/routes.py:137           | user_id              | None                   | Depends(require_role(ADMIN)) | Admin only |
| GET         | `/`                    | root                            | nox_api_v7_fixed.py:128      | None                 | JSON                   | None                         | API root   |
| GET         | `/api/v7/status`       | api_status                      | nox_api_v7_fixed.py:153      | None                 | JSON                   | None                         | Status     |
| GET         | `/api/v7/health`       | health_check                    | nox_api_v7_fixed.py:155      | None                 | JSON                   | None                         | Health     |
| GET         | `/api/v7/metrics/prometheus` | get_prometheus_metrics    | nox_api_v7_fixed.py:165      | None                 | Prometheus text        | None                         | Metrics    |
| GET         | `/health`              | health                          | nox-api/api/nox_api_v23.py:66| request              | JSON                   | Optional Auth                | Health     |
| GET         | `/metrics`             | metrics                         | nox-api/api/nox_api_v23.py:75| current_user         | Prometheus text        | Optional Auth                | Metrics    |
| POST        | `/put`                 | put                             | nox-api/api/nox_api_v23.py:83| path, file           | JSON/status            | Depends(get_current_user)     | Upload     |
| POST        | `/run_py`              | run_py                          | nox-api/api/nox_api_v23.py:92| RunPy, request       | Output                 | Depends(get_current_user)     | Exec py    |

---

## Middleware & Dependencies

### Authentication Middleware
- **get_current_user**: Validates JWT tokens, returns current user
- **require_role**: Role-based access control (ADMIN, USER, etc.)
- **Optional Auth**: Some endpoints allow anonymous access

### Request/Response Patterns
- **Standard Response**: Most endpoints return JSON with consistent structure
- **Error Handling**: Standardized error responses with status codes
- **Pagination**: List endpoints support limit/offset parameters

### File Upload/Processing
- **PUT /put**: File upload with path specification
- **POST /run_py**: Python code execution with request context

---

## API Versioning
- **v7**: Current stable version (`/api/v7/*`)
- **v23**: Development version (`nox-api/api/nox_api_v23.py`)
- **Root**: Unversioned endpoints for health/metrics

---

## Security Features
- JWT-based authentication
- Role-based authorization
- Optional authentication for monitoring endpoints
- File upload restrictions
- Code execution sandboxing

---

**Note**: This inventory is based on available search results and may not include all endpoints. Review actual codebase for complete API documentation.
