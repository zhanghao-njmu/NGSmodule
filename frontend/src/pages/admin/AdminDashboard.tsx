/**
 * Admin Dashboard - System overview and user management
 */
import React, { useEffect, useState } from 'react'
import {
  Card,
  Tag,
  Button,
  Space,
  Modal,
  Form,
  Input,
  Select,
  InputNumber,
  Switch,
  Progress,
  message,
  Popconfirm,
} from 'antd'
import {
  UserOutlined,
  ProjectOutlined,
  ExperimentOutlined,
  ThunderboltOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  DatabaseOutlined,
  EditOutlined,
  DeleteOutlined,
  StopOutlined,
  PlayCircleOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { adminService } from '../../services/admin.service'
import { StatisticCard, DataTable, StatusTag } from '../../components/common'
import type { StatisticItem } from '../../components/common'
import type { User, UserAdminUpdate, SystemStats } from '../../types/admin'
import dayjs from 'dayjs'

const { Option } = Select

export const AdminDashboard: React.FC = () => {
  const [stats, setStats] = useState<SystemStats | null>(null)
  const [users, setUsers] = useState<User[]>([])
  const [loading, setLoading] = useState(false)
  const [editModalOpen, setEditModalOpen] = useState(false)
  const [selectedUser, setSelectedUser] = useState<User | null>(null)
  const [form] = Form.useForm()

  useEffect(() => {
    loadData()
  }, [])

  const loadData = async () => {
    setLoading(true)
    try {
      const [systemStats, userList] = await Promise.all([
        adminService.getSystemStats(),
        adminService.getUsers({ limit: 100 }),
      ])
      setStats(systemStats)
      setUsers(userList)
    } catch (error: any) {
      message.error(`Failed to load data: ${error.message}`)
    } finally {
      setLoading(false)
    }
  }

  const handleEdit = (user: User) => {
    setSelectedUser(user)
    form.setFieldsValue({
      full_name: user.full_name,
      email: user.email,
      organization: user.organization,
      role: user.role,
      is_active: user.is_active,
      storage_quota: Math.round(user.storage_quota / (1024 * 1024 * 1024)), // Convert to GB
    })
    setEditModalOpen(true)
  }

  const handleEditSubmit = async () => {
    if (!selectedUser) return

    try {
      const values = await form.validateFields()
      const updateData: UserAdminUpdate = {
        ...values,
        storage_quota: values.storage_quota * 1024 * 1024 * 1024, // Convert GB to bytes
      }

      await adminService.updateUser(selectedUser.id, updateData)
      message.success('User updated successfully')
      setEditModalOpen(false)
      loadData()
    } catch (error: any) {
      message.error(`Failed to update user: ${error.message}`)
    }
  }

  const handleToggleStatus = async (user: User) => {
    try {
      await adminService.toggleUserStatus(user.id)
      message.success(`User ${user.is_active ? 'deactivated' : 'activated'} successfully`)
      loadData()
    } catch (error: any) {
      message.error(`Failed to toggle user status: ${error.message}`)
    }
  }

  const handleDelete = async (userId: string) => {
    try {
      await adminService.deleteUser(userId)
      message.success('User deleted successfully')
      loadData()
    } catch (error: any) {
      message.error(`Failed to delete user: ${error.message}`)
    }
  }

  const formatBytes = (bytes: number) => {
    const gb = bytes / (1024 * 1024 * 1024)
    return `${gb.toFixed(2)} GB`
  }

  const columns: ColumnsType<User> = [
    {
      title: 'Username',
      dataIndex: 'username',
      key: 'username',
      render: (username) => <span style={{ fontWeight: 500 }}>{username}</span>,
    },
    {
      title: 'Email',
      dataIndex: 'email',
      key: 'email',
    },
    {
      title: 'Role',
      dataIndex: 'role',
      key: 'role',
      width: 100,
      render: (role) => (
        <Tag color={role === 'admin' ? 'red' : 'blue'}>{role.toUpperCase()}</Tag>
      ),
    },
    {
      title: 'Status',
      dataIndex: 'is_active',
      key: 'is_active',
      width: 100,
      render: (is_active) => (
        <StatusTag status={is_active ? 'active' : 'inactive'} />
      ),
    },
    {
      title: 'Storage',
      key: 'storage',
      width: 200,
      render: (_, record) => {
        const percent = (record.storage_used / record.storage_quota) * 100
        return (
          <div>
            <Progress
              percent={Math.round(percent)}
              size="small"
              status={percent > 90 ? 'exception' : 'normal'}
            />
            <div style={{ fontSize: 12, color: '#666' }}>
              {formatBytes(record.storage_used)} / {formatBytes(record.storage_quota)}
            </div>
          </div>
        )
      },
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 120,
      render: (date) => dayjs(date).format('YYYY-MM-DD'),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 180,
      render: (_, record) => (
        <Space>
          <Button
            size="small"
            icon={<EditOutlined />}
            onClick={() => handleEdit(record)}
          >
            Edit
          </Button>
          <Button
            size="small"
            icon={record.is_active ? <StopOutlined /> : <PlayCircleOutlined />}
            onClick={() => handleToggleStatus(record)}
          >
            {record.is_active ? 'Deactivate' : 'Activate'}
          </Button>
          <Popconfirm
            title="Delete user"
            description="Are you sure? This will delete all user data."
            onConfirm={() => handleDelete(record.id)}
            okText="Yes"
            cancelText="No"
          >
            <Button size="small" danger icon={<DeleteOutlined />} />
          </Popconfirm>
        </Space>
      ),
    },
  ]

  const statisticItems: StatisticItem[] = [
    {
      key: 'users',
      title: 'Total Users',
      value: stats?.total_users || 0,
      prefix: <UserOutlined />,
      valueStyle: { color: 'var(--color-primary)' },
      suffix: <div style={{ fontSize: 12, color: 'var(--text-secondary)', marginTop: 8 }}>{stats?.active_users || 0} active</div>,
    },
    {
      key: 'items',
      title: 'Total Projects',
      value: stats?.total_projects || 0,
      prefix: <ProjectOutlined />,
      valueStyle: { color: 'var(--color-success)' },
    },
    {
      key: 'samples',
      title: 'Total Samples',
      value: stats?.total_samples || 0,
      prefix: <ExperimentOutlined />,
      valueStyle: { color: 'var(--color-warning)' },
    },
    {
      key: 'tasks',
      title: 'Total Tasks',
      value: stats?.total_tasks || 0,
      prefix: <ThunderboltOutlined />,
      valueStyle: { color: '#722ed1' },
      suffix: (
        <div style={{ fontSize: 12, marginTop: 8 }}>
          <CheckCircleOutlined style={{ color: 'var(--color-success)' }} /> {stats?.completed_tasks || 0}{' '}
          | <CloseCircleOutlined style={{ color: 'var(--color-error)' }} /> {stats?.failed_tasks || 0}
        </div>
      ),
    },
  ]

  return (
    <div>
      <h2>Admin Dashboard</h2>

      {/* System Statistics */}
      <StatisticCard items={statisticItems} gutter={[16, 16]} />

      {/* Storage Statistics */}
      <Card
        title={
          <>
            <DatabaseOutlined /> Storage Usage
          </>
        }
        style={{ marginBottom: 24 }}
      >
        <Progress
          percent={Math.round(
            ((stats?.total_storage_used || 0) / (stats?.total_storage_quota || 1)) * 100
          )}
          status="active"
        />
        <div style={{ marginTop: 8, fontSize: 14 }}>
          {formatBytes(stats?.total_storage_used || 0)} /{' '}
          {formatBytes(stats?.total_storage_quota || 0)}
        </div>
      </Card>

      {/* User Management Table */}
      <DataTable
        title="User Management"
        columns={columns}
        dataSource={users}
        rowKey="id"
        loading={loading}
        pagination={{ pageSize: 20 }}
        emptyText="No Users"
        emptyDescription="No users have been registered yet"
      />

      {/* Edit User Modal */}
      <Modal
        title="Edit User"
        open={editModalOpen}
        onCancel={() => setEditModalOpen(false)}
        onOk={handleEditSubmit}
        width={600}
      >
        <Form form={form} layout="vertical">
          <Form.Item name="full_name" label="Full Name">
            <Input placeholder="Enter full name" />
          </Form.Item>

          <Form.Item name="email" label="Email" rules={[{ type: 'email' }]}>
            <Input placeholder="Enter email" />
          </Form.Item>

          <Form.Item name="organization" label="Organization">
            <Input placeholder="Enter organization" />
          </Form.Item>

          <Form.Item name="role" label="Role" rules={[{ required: true }]}>
            <Select>
              <Option value="user">User</Option>
              <Option value="admin">Admin</Option>
            </Select>
          </Form.Item>

          <Form.Item name="is_active" label="Active" valuePropName="checked">
            <Switch />
          </Form.Item>

          <Form.Item
            name="storage_quota"
            label="Storage Quota (GB)"
            rules={[{ required: true, type: 'number', min: 1 }]}
          >
            <InputNumber style={{ width: '100%' }} min={1} />
          </Form.Item>
        </Form>
      </Modal>
    </div>
  )
}

export default AdminDashboard
